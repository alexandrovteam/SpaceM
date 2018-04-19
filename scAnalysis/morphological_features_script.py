import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt

#MeasureObjectSizeShape http://cellprofiler.org/manuals/current/MeasureObjectSizeShape.html
#MeasureObjectIntensity http://cellprofiler.org/manuals/current/MeasureObjectIntensity.html
#MeasureObjectNeighbors http://cellprofiler.org/manuals/current/MeasureObjectNeighbors.html
features_OI = ['AreaShape_Area',
               'AreaShape_Compactness',
               'AreaShape_Eccentricity',
               'AreaShape_EulerNumber',
               'Intensity_IntegratedIntensity_OrigRed_highDR',
               'Intensity_MADIntensity_OrigRed_highDR',
               'Intensity_MaxIntensity_OrigRed_highDR',
               'Intensity_MeanIntensity_OrigRed_highDR',
               'Intensity_MedianIntensity_OrigRed_highDR',
               'Intensity_MinIntensity_OrigRed_highDR',
               'Intensity_StdIntensity_OrigRed_highDR',
               'Location_Center_X',
               'Location_Center_Y',
               'Neighbors_FirstClosestDistance_Adjacent',
               'Neighbors_NumberOfNeighbors_Adjacent',
               'Neighbors_PercentTouching_Adjacent',
               'Neighbors_SecondClosestDistance_Adjacent']

features_OI_names = ['Area',
               'Compactness',
               'Eccentricity',
               'EulerNumber',
               'Intensity SUM',
               'Intensity MAD',
               'Intensity MAX',
               'Intensity MEAN',
               'Intensity MEDIAN',
               'Intensity_MIN',
               'Intensity_STD',
               'Location X',
               'Location Y',
               'FirstClosestDistance',
               'NumberOfNeighbors',
               'PercentTouching',
               'SecondClosestDistance']

MFA = 'D:/Experiments/Hepa_June/Data/transformation_1/TNFa_2.2/Analysis/'
CD = 0.75
CD_fname = MFA + 'tSNE/CD={}_tSNEdata/'.format(CD)
tsne_input, fluoCell_input, fluoMarks_input, keys_input, tsne_sf, tsne_adduct, tsne_data, best_corrS, \
CD, pmi, mbi, tsne_input_norm, fluoMarks_input_norm, overlapCellMarks_input, cell_area_input \
    = np.load(CD_fname + 'tsne_inputs_CD=' + str(CD) + '.npy')

feat_df = pd.read_csv("D:\Experiments\Hepa_June\Data/transformation_1\TNFa_2.2\Analysis\CellProfilerAnalysis\Cells.csv")
mol_df = pd.DataFrame(data=tsne_input, columns=tsne_sf)
ObjN = np.array([int(i) for i in keys_input])
plt.scatter(cell_area_input, feat_df.iloc[ObjN-1,np.ones(ObjN.shape)*2].as_matrix()[:,0])
plt.xlabel('Cell area calculated from custom script', fontsize=15)
plt.ylabel('Cell area calculated from CellProfiler', fontsize=15)
feat_df_selected = feat_df[features_OI]
feat_df_selected.columns = features_OI_names
df = pd.concat([feat_df_selected.iloc[ObjN-1,:].reset_index(), mol_df.reset_index()], axis=1)
# corr = df.corr(method='spearman')
df.corr(method='spearman').to_csv("D:\Experiments\Hepa_June\Data/transformation_1\TNFa_2.2\Analysis\CellProfilerAnalysis/mol_morph_features_ClusterGrammer_input.txt", sep="\t")

plt.figure()
sns.regplot(df.NumberOfNeighbors,df.C19H24N2O2)

plt.figure()
sns.regplot(df.NumberOfNeighbors,df.C44H74O13P2)

sns.pairplot(df,
             x_vars=["NumberOfNeighbors"],
             y_vars=["C19H24N2O2", "C44H74O13P2"],
             size=5, aspect=.8,
             kind="reg")
common_mol_sf = ['C18H32O2', 'C20H32O2', 'C21H37O6P', 'C21H39O6P', 'C23H46NO7P', 'C25H44NO7P', 'C37H69O8P', 'C39H69O8P', 'C39H71O8P', 'C39H73O8P', 'C39H76NO8P', 'C41H71O8P', 'C41H73O8P', 'C41H74NO8P', 'C41H78NO8P','C43H81O13P', 'C45H66O5', 'C45H68O5', 'C45H86O13P', 'C46H79NO10P',
          'C47H70O5', 'C48H81NO10P', 'C48H83NO10P', 'C55H72N4O6']
common_mol_name = ['Linoleic acid', 'Arachidonic acid','CPA(18:2)','CPA(18:1)','LysoPE(18:1)','LysoPE(20:4)','PA(34:2)', 'PA(36:4)','PA(36:3)','PA(36:2)','PE(34:1)','PA(38:5)','PA(38:4)','PE(36:4)','PE(36:2)','PI(34:1)',
              'DG(42:11)','DG(42:10)','PI(36:0)','PS(40:5)','DG(44:11)','PS(42:6)','PS(42:5)','pheophytin b']
for i,sf in enumerate(common_mol_sf):
    sns.jointplot(x="NumberOfNeighbors", y=sf, data=df, kind="reg")
    plt.savefig('D:\Experiments\Hepa_June\Data/transformation_1\TNFa_2.2\Analysis/'
                'CellProfilerAnalysis/neighbor_mol_corrplots_common_mol/' + common_mol_name[i].replace(':', '_') + '.png', dpi=200)
    plt.close('all')
    print(i)