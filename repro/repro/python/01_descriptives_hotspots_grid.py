"""Descriptives, beta-binomial river/season tables, logistic regression, Moran's I, Getis-Ord Gi* hotspots, and the Kaggle point-grid check.
Run from the repository root:  python python/01_descriptives_hotspots_grid.py
Inputs: data/Water quality.csv and data/Pollution_Intensity_Coordinates_Updated.csv (override with HYPOXIA_DATA_FILE and HYPOXIA_GRID_FILE). Outputs go to results/ (stations.csv is needed by 02)."""
import pandas as pd, numpy as np, json, sys
from scipy.stats import beta, spearmanr, kruskal, mannwhitneyu
from scipy.spatial import cKDTree
import statsmodels.formula.api as smf
import os
DATA_FILE=os.environ.get('HYPOXIA_DATA_FILE','data/Water quality.csv')   # Sinha et al. (2026), Mendeley Data
GRID_FILE=os.environ.get('HYPOXIA_GRID_FILE','data/Pollution_Intensity_Coordinates_Updated.csv')   # Kaggle Bangladesh-River-Pollution-GIS point grid
os.makedirs('results',exist_ok=True)
w=pd.read_csv(DATA_FILE,encoding='utf-8-sig'); w.columns=[c.replace('\n',' ') for c in w.columns]
w=w.rename(columns={'DO (mg/L)':'DO','COD (mg/L)':'COD','Nitrate (mg/L)':'NO3','Phosphate (mg/L)':'PO4','Turbidity (NTU)':'Turb','Temp (°C)':'Temp','Rainfall (mm)':'Rain','River/Stream':'River'})
os.chdir('results')
w['River']=w.River.str.strip().replace({'Tongi khal':'Tongi Khal'})
w['Season']=w.Season.str.replace('-',' ').str.strip().replace({'Pre monsoon':'Pre-monsoon','Post monsoon':'Post-monsoon'})
w['date']=pd.to_datetime(w.Date); w['year']=w.date.dt.year
w['st']=w.River+' | '+w.Location
R={}
R['n']=len(w); R['stations']=w.st.nunique(); R['rivers']=w.River.nunique()
R['cod70_share']=round((w.COD==70).mean(),3)
R['dup_rows']=int(w.duplicated().sum())
# thresholds ECR 2023 fisheries
w['lt5']=(w.DO<5).astype(int); w['lt2']=(w.DO<2).astype(int)
w['fail_cod']=w.COD>50; w['fail_ph']=(w.pH<6)|(w.pH>9)
w['habit_fail']=(w.lt5==1)|w.fail_cod|w.fail_ph
R['p_lt5']=round(w.lt5.mean(),3); R['p_lt2']=round(w.lt2.mean(),3); R['p_codfail']=round(w.fail_cod.mean(),3); R['p_habitfail']=round(w.habit_fail.mean(),3)
R['DO_median']=float(w.DO.median())
def post(k,n): a,b=k+0.5,n-k+0.5; return a/(a+b), beta.ppf(0.025,a,b), beta.ppf(0.975,a,b)
def tab(by,col='lt2'):
    g=w.groupby(by)[col].agg(['sum','size']).reset_index()
    g[['p','lo','hi']]=g.apply(lambda r: pd.Series(post(r['sum'],r['size'])),axis=1)
    return g
rs=tab(['River','Season']); rs.to_csv('river_season_lt2.csv',index=False)
r5=tab(['River'],'lt5'); r2=tab(['River'],'lt2')
riv=r5.merge(r2,on='River',suffixes=('_5','_2'))
riv['DOmed']=riv.River.map(w.groupby('River').DO.median())
riv=riv.sort_values('p_2',ascending=False); riv.to_csv('river.csv',index=False)
se=tab(['Season'],'lt2'); se5=tab(['Season'],'lt5'); 
R['season_lt2']=se.set_index('Season')[['p','lo','hi','size']].round(3).to_dict('index')
R['season_lt5']=se5.set_index('Season')[['p','size']].round(3).to_dict('index')
R['year_lt2']=w.groupby('year').lt2.mean().round(3).to_dict()
R['year_DOmed']=w.groupby('year').DO.median().round(2).to_dict()
R['kw_season_DO']=kruskal(*[g.DO for _,g in w.groupby('Season')]).pvalue
# logistic regression for P(DO<2)
m=smf.logit('lt2 ~ C(Season,Treatment("Winter")) + C(River,Treatment("Shitalakhya")) + Temp + np.log1p(Rain) + year',data=w[w.River.map(w.River.value_counts())>=40]).fit(disp=0)
R['logit_pseudoR2']=round(m.prsquared,3)
co=pd.DataFrame({'OR':np.exp(m.params),'lo':np.exp(m.conf_int()[0]),'hi':np.exp(m.conf_int()[1]),'p':m.pvalues}).round(3)
co.to_csv('logit_lt2.csv')
# Monte Carlo: station-level stochastic risk via bootstrap of seasonal draws (prob that a random monthly sample in a season is hypoxic, with uncertainty)
st=w.groupby('st').agg(River=('River','first'),Location=('Location','first'),lat=('Latitude','mean'),lon=('Longitude','mean'),n=('DO','size'),k2=('lt2','sum'),k5=('lt5','sum'),DOmed=('DO','median'),COD=('COD','median'),WQI=('WQI','mean')).reset_index()
st[['p2','lo2','hi2']]=st.apply(lambda r: pd.Series(post(r.k2,r.n)),axis=1)
st['p5']=(st.k5+0.5)/(st.n+1)
# spatial autocorrelation Moran's I on station p2 (inverse distance, 5 km band)
la0=23.8; xy=np.c_[st.lat*111.0, st.lon*111.0*np.cos(np.radians(la0))]
D=np.sqrt(((xy[:,None,:]-xy[None,:,:])**2).sum(-1)); W=np.where((D>0)&(D<=5),1/np.maximum(D,0.1),0); 
def moran(x,W):
    z=x-x.mean(); s0=W.sum(); I=len(x)/s0*(z@W@z)/(z@z); 
    rng=np.random.default_rng(1); sims=[]
    for _ in range(999):
        zp=rng.permutation(z); sims.append(len(x)/s0*(zp@W@zp)/(zp@zp))
    return I,(np.sum(np.array(sims)>=I)+1)/1000
st_use=st[st.n>=5].reset_index(drop=True)
xy2=np.c_[st_use.lat*111.0, st_use.lon*111.0*np.cos(np.radians(la0))]; D2=np.sqrt(((xy2[:,None,:]-xy2[None,:,:])**2).sum(-1)); W2=np.where((D2>0)&(D2<=5),1/np.maximum(D2,0.1),0)
I,p=moran(st_use.p2.values,W2); R['moran_p2']=(round(I,3),round(p,3)); R['moran_n']=len(st_use); R['moran_isolated']=int((W2.sum(1)==0).sum())
# Getis-Ord Gi* hotspots
Wg=W2.copy(); np.fill_diagonal(Wg,1)
x=st_use.p2.values; n=len(x); xb=x.mean(); S=np.sqrt((x**2).mean()-xb**2)
wi=Wg.sum(1); wi2=(Wg**2).sum(1)
Gi=(Wg@x - xb*wi)/(S*np.sqrt((n*wi2-wi**2)/(n-1)))
st_use['Gi']=Gi
# grid validation
g=pd.read_csv(GRID_FILE,encoding='utf-8-sig')
gxy=np.c_[g.Latitude*111.0, g.Longitude*111.0*np.cos(np.radians(la0))]; t=cKDTree(gxy)
d,i=t.query(np.c_[st.lat*111.0, st.lon*111.0*np.cos(np.radians(la0))]); st['grid_near']=g.Pollution_Intensity_Value.values[i]; st['grid_d']=d
V={}
for rad in (0.5,1.0,1.5):
    vals=[g.Pollution_Intensity_Value.values[j].mean() if len(j) else np.nan for j in t.query_ball_point(np.c_[st.lat*111.0, st.lon*111.0*np.cos(np.radians(la0))],rad)]
    st[f'grid_{rad}']=vals
    s=st.dropna(subset=[f'grid_{rad}'])
    V[rad]={m:[round(v,3) for v in spearmanr(s[f'grid_{rad}'],s[m])] for m in ('p2','DOmed','COD','WQI')}; V[rad]['n']=len(s)
R['grid_valid']=V
# leave-one-river-out stability of near-point rho with p2
lo=[]
for rv in st.River.unique():
    s=st[st.River!=rv]; lo.append(spearmanr(s.grid_near,s.p2)[0])
R['grid_lorO_range']=(round(min(lo),3),round(max(lo),3))
# grid internals
R['grid_class_counts']=g.Pollution_Intensity_Class.value_counts().to_dict()
gs=g.sort_values(['Latitude','Longitude']); R['grid_spacing_deg']=float(np.round(np.median(np.diff(np.sort(g.Longitude.unique()))),5))
gx=np.c_[g.Latitude*111.0, g.Longitude*111.0*np.cos(np.radians(la0))]; tg=cKDTree(gx); nn=tg.query(gx,k=9)
R['grid_nn_rho']=round(spearmanr(g.Pollution_Intensity_Value, g.Pollution_Intensity_Value.values[nn[1][:,1:]].mean(1))[0],3)
st.to_csv('stations.csv',index=False); st_use.to_csv('stations_gi.csv',index=False)
R['hotspots']=st_use[st_use.Gi>1.96][['River','Location','n','p2','Gi']].round(2).to_dict('records')
R['coldspots']=st_use[st_use.Gi<-1.96][['River','Location','n','p2','Gi']].round(2).to_dict('records')
R['top_stations']=st[st.n>=10].sort_values('p2',ascending=False).head(8)[['River','Location','n','k2','p2','lo2','hi2','DOmed']].round(2).to_dict('records')
json.dump(R,open('results.json','w'),indent=1,default=float)
print(json.dumps(R,indent=1,default=float)); print(riv.round(3).to_string()); print(co.to_string())
