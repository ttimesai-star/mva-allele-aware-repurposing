"""Sequence- and structure-based effect estimates for BUBR1 missense variants (public data only).

Inputs are public: UniProt O60566 sequence, the AlphaFold DB model AF-O60566-F1, ClinVar BUB1B missense sets
(benign/likely benign) and the published MVA1 missense alleles. The proband's variant enters only as its protein
change (N1002K). Heavy (CPU ~30 min for both ESM models; ~2 GB RAM per model):

    python src/variant_effect_models.py esm --workdir D:/work         # ESM-1v + ESM-2 masked-marginal LLR, AlphaFold context
    python src/variant_effect_models.py ddg --workdir D:/work --thermompnn D:/ThermoMPNN   # ThermoMPNN ddG, kinase-domain crop

Writes data/processed/esm_scores.json and data/processed/thermompnn_kd.csv (copied from --workdir).
Models: ESM-1v (Meier et al. 2021, model 1 of 5), ESM-2 650M (Lin 2023), ThermoMPNN (Dieckhaus 2024; positive ddG =
destabilising in this output), AlphaFold DB v6 model (Jumper 2021; Varadi 2024).
"""
import sys
MODE = sys.argv[1] if len(sys.argv) > 1 else "esm"
if MODE == "esm":
    # Public-data only: UniProt O60566 sequence + AlphaFold model; scores named missense substitutions.
    import json, requests, os, math, numpy as np, torch, time
    from transformers import AutoTokenizer, EsmForMaskedLM
    W=sys.argv[sys.argv.index('--workdir')+1]; os.makedirs(W,exist_ok=True)
    seq=''.join(requests.get('https://rest.uniprot.org/uniprotkb/O60566.fasta',timeout=60).text.split('\n')[1:])
    print('len',len(seq), 'res1002', seq[1001], 'res737', seq[736])
    three={'Ala':'A','Arg':'R','Asn':'N','Asp':'D','Cys':'C','Gln':'Q','Glu':'E','Gly':'G','His':'H','Ile':'I','Leu':'L','Lys':'K','Met':'M','Phe':'F','Pro':'P','Ser':'S','Thr':'T','Trp':'W','Tyr':'Y','Val':'V'}
    sets=json.load(open(sys.argv[sys.argv.index('--sets')+1] if '--sets' in sys.argv else 'data/processed/clinvar_missense_sets.json'))
    def p(v):
        if v[1].islower(): return three[v[:3]],int(v[3:-3]),three[v[-3:]]
        return v[0],int(v[1:-1]),v[-1]
    benign=[p(v) for v in sets['B']]
    lit_low=[p(v) for v in ['R727C','R814H','L844F','I909T','L1012P']]
    lit_norm=[p(v) for v in ['Y155C','Q921H']]
    query=[('N',1002,'K')]
    AA='ACDEFGHIKLMNPQRSTVWY'
    def window(pos):
        L=len(seq); 
        if L<=1022: return 0
        s=max(0,min(pos-511,L-1022)); return s
    res={}
    for mname in ['facebook/esm1v_t33_650M_UR90S_1','facebook/esm2_t33_650M_UR50D']:
        t0=time.time()
        tok=AutoTokenizer.from_pretrained(mname); model=EsmForMaskedLM.from_pretrained(mname).eval()
        cache={}
        def llr(pos):
            if pos in cache: return cache[pos]
            s=window(pos); sub=seq[s:s+1022]; i=pos-1-s
            ids=tok(sub,return_tensors='pt')['input_ids']
            ids[0,i+1]=tok.mask_token_id
            with torch.no_grad(): lg=model(ids).logits[0,i+1]
            lp=torch.log_softmax(lg,-1)
            d={a:lp[tok.convert_tokens_to_ids(a)].item() for a in AA}
            cache[pos]=d; return d
        def score(w,pos,m):
            assert seq[pos-1]==w,(w,pos,seq[pos-1])
            d=llr(pos); return d[m]-d[w]
        r={'query':{f'{w}{p}{m}':score(w,p,m) for w,p,m in query},
           'lit_low':{f'{w}{p}{m}':score(w,p,m) for w,p,m in lit_low},
           'lit_normal':{f'{w}{p}{m}':score(w,p,m) for w,p,m in lit_norm},
           'clinvar_benign':{f'{w}{p}{m}':score(w,p,m) for w,p,m in benign}}
        # saturation over kinase domain 766-1050 for percentile
        allkd=[]
        for pos in range(766,1051):
            d=llr(pos); w=seq[pos-1]
            allkd+= [d[a]-d[w] for a in AA if a!=w]
        q=list(r['query'].values())[0]
        r['kd_percentile_of_query']=float(np.mean(np.array(allkd)>=q))
        r['site1002_all']={a:llr(1002)[a]-llr(1002)['N'] for a in AA}
        res[mname]=r; print(mname,'done',round(time.time()-t0),'s', json.dumps(r['query']), 'pct',r['kd_percentile_of_query'])
    json.dump(res,open(W+'/esm_scores.json','w'),indent=1)
    print(json.dumps(res))
    # AlphaFold structure context
    af=requests.get('https://alphafold.ebi.ac.uk/api/prediction/O60566',timeout=60).json()[0]
    pdb=requests.get(af['pdbUrl'],timeout=120).text; open(W+'/AF-O60566.pdb','w').write(pdb)
    print('AF version',af.get('latestVersion'),af['pdbUrl'])
    from Bio.PDB import PDBParser, ShrakeRupley, NeighborSearch
    st=PDBParser(QUIET=True).get_structure('a',W+'/AF-O60566.pdb')
    sr=ShrakeRupley(); sr.compute(st[0],level='R')
    maxasa={'A':129,'R':274,'N':195,'D':193,'C':167,'Q':225,'E':223,'G':104,'H':224,'I':197,'L':201,'K':236,'M':224,'F':240,'P':159,'S':155,'T':172,'W':285,'Y':263,'V':174}
    ch=st[0]['A']; atoms=[a for a in st.get_atoms()]; ns=NeighborSearch(atoms)
    from Bio.PDB.Polypeptide import three_to_index, index_to_one
    def ctx(pos):
        r=ch[pos]; one=index_to_one(three_to_index(r.get_resname()))
        rsa=r.sasa/maxasa[one]; pl=r['CA'].get_bfactor()
        nb=sorted({a.get_parent().id[1] for a in ns.search(r['CA'].coord,8.0)}-{pos})
        hb=[]
        for at in r:
            if at.element in('N','O') and at.get_id() not in('N','O'):
                for b in ns.search(at.coord,3.5):
                    if b.element in('N','O') and b.get_parent().id[1]!=pos: hb.append((at.get_id(),b.get_parent().get_resname()+str(b.get_parent().id[1]),b.get_id(),round(float(np.linalg.norm(at.coord-b.coord)),2)))
        return dict(pos=pos,res=one,pLDDT=pl,RSA=round(rsa,3),n_CA_neighbours_8A=len(nb),sidechain_polar_contacts=hb)
    for pos in [1002,1012,909,844,814,727]:
        print('CTX',json.dumps(ctx(pos)))
else:
    # ThermoMPNN ddG on the public AlphaFold model of BUBR1 (O60566), kinase-domain crop. Public data only.
    import subprocess, os, types
    W=sys.argv[sys.argv.index('--workdir')+1]
    subprocess.run([sys.executable,'-m','pip','install','-q','torchmetrics'],capture_output=True)
    sys.modules['wandb']=types.ModuleType('wandb')
    import torch
    _tl=torch.load
    def tl(*a,**k):
        k['weights_only']=False; k.setdefault('map_location','cpu'); return _tl(*a,**k)
    torch.load=tl
    TD=sys.argv[sys.argv.index('--thermompnn')+1]
    sys.path.insert(0,TD); sys.path.insert(0,TD+'/analysis')
    from omegaconf import OmegaConf
    import pandas as pd
    # crop AF model to residues 720-1050 (kinase domain + N-terminal extension)
    src=W+'/AF-O60566.pdb'; dst=W+'/BUBR1_KD_720_1050.pdb'
    with open(src) as f, open(dst,'w') as g:
        for l in f:
            if l.startswith('ATOM') and 720<=int(l[22:26])<=1050: g.write(l)
        g.write('END\n')
    cfg=OmegaConf.load(TD+'/local.yaml')
    cfg=OmegaConf.merge(cfg,OmegaConf.create({'training':{'num_workers':1,'learn_rate':0.001,'epochs':100,'lr_schedule':True},
      'model':{'hidden_dims':[64,32],'subtract_mut':True,'num_final_layers':2,'freeze_weights':True,'load_pretrained':True,'lightattn':True,'lr_schedule':True}}))
    cfg.platform.thermompnn_dir=TD; cfg.platform.accel='cpu'
    from train_thermompnn import TransferModelPL
    from protein_mpnn_utils import alt_parse_PDB
    from datasets import Mutation
    from SSM import get_ssm_mutations
    model=TransferModelPL.load_from_checkpoint(TD+'/models/thermoMPNN_default.pt',cfg=cfg,map_location='cpu').model.eval()
    pdb=alt_parse_PDB(dst,'A')
    muts=get_ssm_mutations(pdb[0])
    ml=[]
    for m in muts:
        if m is None: ml.append(None); continue
        m=m.strip(); ml.append(Mutation(position=int(m[1:-1]),wildtype=m[0],mutation=m[-1],ddG=None,pdb=pdb[0]['name']))
    with torch.no_grad():
        pred,_=model(pdb,ml)
    rows=[]
    for mu,o in zip(ml,pred):
        if mu is not None: rows.append((mu.position,mu.wildtype,mu.mutation,o['ddG'].cpu().item()))
    df=pd.DataFrame(rows,columns=['idx','wt','mut','ddG'])
    # ThermoMPNN positions are 0-based indices into the parsed chain sequence -> map to UniProt numbering
    seq=pdb[0]['seq_chain_A'] if 'seq_chain_A' in pdb[0] else pdb[0]['seq']
    print('parsed len',len(seq),'first res',seq[:10])
    off=720 - (df.idx.min())
    df['pos']=df.idx+off
    df.to_csv(W+'/thermompnn_kd.csv',index=False)
    def get(p,w,m):
        r=df[(df.pos==p)&(df.mut==m)]; 
        return None if r.empty else (r.wt.iloc[0], round(float(r.ddG.iloc[0]),3))
    print('README sign lines:',[l for l in open(TD+'/README.md',encoding='utf8') if 'stabiliz' in l.lower()][:5])
    for v in [('N',1002,'K'),('L',1012,'P'),('I',909,'T'),('L',844,'F'),('R',814,'H'),('R',727,'C'),('Q',921,'H'),('A',784,'T'),('D',818,'E'),('F',822,'L'),('T',859,'A')]:
        print(v, get(v[1],v[0],v[2]))
    import numpy as np
    allv=df[df.wt!=df.mut].ddG.values
    n=get(1002,'N','K')[1]
    print('fraction of KD substitutions more destabilising than N1002K:', float(np.mean(allv<n)), 'median all', float(np.median(allv)))
    print('site 1002 all:', df[df.pos==1002][['mut','ddG']].round(2).to_dict('records'))
