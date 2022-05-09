from aminoAcids import *
from main import *
import matplotlib.pyplot as plt
import time
import numpy as np
import pandas as pd
import prody as prdy
import seaborn as sns

def scale(ag, fname):
    coords = ag.getCoords()
    coords *= 3.8
    # coords = np.multiply(coords, 3.8)
    ag.setCoords(coords)
    prdy.writePDB(fname,ag)

    return ag

def CompareModels(hp_mc, hp_sa, hp_mc_scaled, hp_sa_scaled):
    # m1 = prdy.parsePDB(f_Quark)
    # m1 = m1.select('calpha')
    # prdy.writePDB('model1-calpha', m1)    # save calpha file for VMD
    quark = prdy.parsePDB(f'{f_QUARK}')
    quark = quark.calpha
    # print(f'quark.numAtoms()={quark.numAtoms()}; should be 141')

    # data = prdy.parsePDB('1si4.pdb')
    # data = data.select('calpha')
    # prdy.writePDB('data-calpha', data)    # save calpha file for VMD
    data = prdy.parsePDB(f'{f_1si4}')
    data = data.calpha
    # print(f'data.numAtoms()={data.numAtoms()}; should be 574')

    trans_quark, quark_chA, data_chA, sequid, overlap = prdy.matchAlign(quark, data)
    RMSD_score_1 = prdy.calcRMSD(quark_chA, data_chA)
    print(f'RMSD scores corresponding to prdy.calcRMSD(quark_chA, data_chA)={RMSD_score_1}')

    trans_hp_mc, hp_mc_chA, data_chA, sequid, overlap = prdy.matchAlign(hp_mc, data)
    RMSD_score_2 = prdy.calcRMSD(hp_mc_chA, data_chA)
    print(f'RMSD scores corresponding to prdy.calcRMSD(hp_mc_chA, data_chA)={RMSD_score_2}')

    trans_hp_mc_scaled, hp_mc_scaled_chA, data_chA, sequid, overlap = prdy.matchAlign(hp_mc_scaled, data)
    RMSD_score_3 = prdy.calcRMSD(hp_mc_scaled_chA, data_chA)
    print(f'RMSD scores corresponding to prdy.calcRMSD(hp_sa_scaled_chA, data_chA)={RMSD_score_3}')

    trans_hp_sa, hp_sa_chA, data_chA, sequid, overlap = prdy.matchAlign(hp_sa, data)
    RMSD_score_4 = prdy.calcRMSD(hp_sa_chA, data_chA)
    print(f'RMSD scores corresponding to prdy.calcRMSD(hp_sa_chA, data_chA)={RMSD_score_4}')

    trans_hp_sa_scaled, hp_sa_scaled_chA, data_chA, sequid, overlap = prdy.matchAlign(hp_sa_scaled, data)
    RMSD_score_5 = prdy.calcRMSD(hp_sa_scaled_chA, data_chA)
    print(f'RMSD scores corresponding to prdy.calcRMSD(hp_sa_scaled_chA, data_chA)={RMSD_score_5}')


    with open(f'RMSD_quark-mc-sa({numstarts}trials)({save_num}).txt', 'w') as fout:
        fout.write(f'RMSD scores corresponding to prdy.calcRMSD(quark_chA, data_chA)={RMSD_score_1}\n')
        fout.write(f'RMSD scores corresponding to prdy.calcRMSD(hp_mc_chA, data_chA)={RMSD_score_2}\n')
        fout.write(f'RMSD scores corresponding to prdy.calcRMSD(hp_mc_scaled_chA, data_chA)={RMSD_score_3}\n')
        fout.write(f'RMSD scores corresponding to prdy.calcRMSD(hp_sa_chA, data_chA)={RMSD_score_4}\n')
        fout.write(f'RMSD scores corresponding to prdy.calcRMSD(hp_sa_scaled_chA, data_chA)={RMSD_score_5}\n')
    
def CompareScoring(energies_sa, allFolds):
    # output all SA folds as pdb files
    rmsd_scores = []
    x = []
    
    data = prdy.parsePDB(f'{f_1si4}')
    data = data.calpha
    
    for i in range(len(allFolds)):
        output_pdb(allFolds[i],f'All_SA_Folds(i)')
        curr_fold = prdy.parsePDB('All_SA_Folds(i).pdb')
        
        trans_hp, hp_chA, data_chA, sequid, overlap = prdy.matchAlign(curr_fold, data)
        rmsd = prdy.calcRMSD(hp_chA, data_chA)
        rmsd_scores.append(rmsd)
        x.append(i+1)
    
    corr = np.corrcoef(energies_sa, rmsd_scores)
    print(f'correlation coefficient =', corr[0, 1])
    
    with open(f'energy_RMSD-coorcoeff({numstarts}trials)({save_num}).txt', 'w') as fout:
        fout.write(f'Correlation Coefficient\n')
        fout.write(f'np.corrcoef(energies_sa, rmsd_scores) = {corr[0,1]}\n')
    
    fig = plt.figure()
    df = pd.DataFrame({'energy': energies_sa, 'RMSD': rmsd_scores})
    
    g = sns.JointGrid(data=df, x='energy', y='RMSD')
    g = g.plot_joint(sns.regplot, color="xkcd:muted blue")
    g = g.plot_marginals(sns.distplot, color="xkcd:bluey grey")
    g.ax_joint.text(-40, 16, f'r = {corr[0,1]}', fontstyle='italic')
    plt.tight_layout()
    plt.savefig(f'correlation-regression({numstarts}trials)({save_num})')


def output_distribution(x, fname):
    fig = plt.figure(figsize=(10,7))
    sns.distplot(x)
    plt.savefig(f'{fname}')

def output_distributions(energies_mc, energies_sa):
    fig = plt.figure(figsize=(10,7))
    df = pd.DataFrame({'energies_mc': energies_mc, 'energies_sa':energies_sa})
    sns.histplot(df, element='step')
    plt.savefig(f'distributions(141aa{numstarts}trials)({save_num})')

def output_pdb(aaSeq, fname):
    
    # template for .pdb output required
    bb = prdy.parsePDB('templatePDB.pdb')
    ag = bb.select('calpha')

    coords = []
    for aa in aaSeq:
        coords += [[aa.r, aa.c, aa.z]]
    ag.setCoords(np.array(coords))

    prdy.writePDB(fname,ag)
    
    return ag

def output_graph(x, y, x_label, y_label, title, fname):
    fig = plt.figure(figsize =(10,7))
    plt.scatter(x, y)
    plt.xlabel(f'{x_label}')
    plt.ylabel(f'{y_label}')
    plt.title(f"{title}")
    plt.savefig(f'{fname}')

def output_plot(x, y, x_label, y_label, title, fname):
    fig = plt.figure(figsize =(10,7))
    plt.plot(x, y)
    plt.xlabel(f'{x_label}')
    plt.ylabel(f'{y_label}')
    plt.title(f"{title}")
    plt.savefig(f'{fname}')

def output_overlay_graph(x, all_y, x_label, y_label, title, fname):
    fig = plt.figure(figsize =(10,7))
    for yi in all_y:
        # print(f'len(i)= {len(i)}; len(x) = {len(x)}')
        if len(yi) < len(x):
            tail = [yi[-1]]*(len(x)-len(yi))
            yi= np.append(yi, [tail])
        # print(f'len(i)= {len(i)}; len(x) = {len(x)}')

        plt.plot(x, yi)
    plt.xlabel(f'{x_label}')
    plt.ylabel(f'{y_label}')
    plt.title(f"{title}")
    plt.savefig(f'{fname}')

def Plot_ProteinLength_vs_Runtime_MC(protein):
    x = []
    start = time.time()
    protein20 = protein[:20]
    aaSeq20 = foldRandom(protein20)
	# for i in range(len(aaSeq20)):
	# 	print(aaSeq20[i].type, aaSeq20[i].r, aaSeq20[i].c)
    energy_20 = energy(aaSeq20, 1)
    print(f'energy for this arrangment of 20aa={energy_20}')
    t = (time.time()-start)/60
    print(f'20aa took {t} min to run')
    x.append(t)

    
    start = time.time()
    protein30 = protein[:30]
    aaSeq30 = foldRandom(protein30)
    energy_30 = energy(aaSeq30, 1)
    print(f'energy for this arrangment of 30aa={energy_30}')
    t = (time.time()-start)/60
    print(f'30aa took {t} min to run')
    x.append(t)

    
    start = time.time()
    protein50 = protein[:50]
    aaSeq50 = foldRandom(protein50)
    energy_50 = energy(aaSeq50, 1)
    print(f'energy for this arrangment of 50aa={energy_50}')
    t = (time.time()-start)/60
    print(f'50aa took {t} min to run')
    x.append(t)
    
    start = time.time()
    protein80 = protein[:80]
    aaSeq80 = foldRandom(protein80)
    energy_80 = energy(aaSeq80, 1)
    print(f'energy for this arrangment of 80aa={energy_80}')
    t = (time.time()-start)/60
    print(f'80aa took {t} min to run')
    x.append(t)
    
    start = time.time()
    aaSeq = foldRandom(protein)
    energy_141 = energy(aaSeq, 1)
    print(f'energy for this arrangment of 141aa={energy_141}')
    t = (time.time()-start)/60
    print(f'{len(aaSeq)}aa took {t} min to run')
    x.append(t)

    # with open('CoV2SpikeProteinSeq.txt') as f:
    #     fin = f.read()
    #     seq = fin.strip().split('\n')
	# # print(seq[0])
    # CoV2SpikeProtein = HP_chain(amino_acid_dict,seq[0],True)
	# # print(CoV2SpikeProtein)
    
    # start = time.time()
    # aaSeqCoV = foldRandom(CoV2SpikeProtein)
    # t = (time.time()-start)
    # print(f'{len(aaSeqCoV)}aa took {(time.time()-start)} seconds to run')
    # x.append(t)

    # i = [20, 30, 50, 80, 141, 1281]
    i = [20, 30, 50, 80, 141]
    output_graph(i, x, "len(protein)","time(sec)","Time vs Protein size (including energy calculation)", f'time-vs-protein-size-MC({save_num})')

def Plot_ProteinLength_vs_Runtime_SA(protein):
    x = []
    start = time.time()
    protein20 = protein[:20]
    aaSeq20 = SA_Sampler(protein20)
	# for i in range(len(aaSeq20)):
	# 	print(aaSeq20[i].type, aaSeq20[i].r, aaSeq20[i].c)
    energy_20 = energy(aaSeq20, 1)
    print(f'energy for this arrangment of 20aa={energy_20}')
    t = (time.time()-start)/60
    print(f'20aa took {t} min to run')
    x.append(t)

    
    start = time.time()
    protein30 = protein[:30]
    aaSeq30 = SA_Sampler(protein30)
    energy_30 = energy(aaSeq30, 1)
    print(f'energy for this arrangment of 30aa={energy_30}')
    t = (time.time()-start)/60
    print(f'30aa took {t} min to run')
    x.append(t)

    
    start = time.time()
    protein50 = protein[:50]
    aaSeq50 = SA_Sampler(protein50)
    energy_50 = energy(aaSeq50, 1)
    print(f'energy for this arrangment of 50aa={energy_50}')
    t = (time.time()-start)/60
    print(f'50aa took {t} min to run')
    x.append(t)
    
    start = time.time()
    protein80 = protein[:80]
    aaSeq80 = SA_Sampler(protein80)
    energy_80 = energy(aaSeq80, 1)
    print(f'energy for this arrangment of 80aa={energy_80}')
    t = (time.time()-start)/60
    print(f'80aa took {t} min to run')
    x.append(t)
    
    start = time.time()
    aaSeq = SA_Sampler(protein)
    energy_141 = energy(aaSeq, 1)
    print(f'energy for this arrangment of 141aa={energy_141}')
    t = (time.time()-start)/60
    print(f'{len(aaSeq)}aa took {t} min to run')
    x.append(t)

    # with open('CoV2SpikeProteinSeq.txt') as f:
    #     fin = f.read()
    #     seq = fin.strip().split('\n')
	# # print(seq[0])
    # CoV2SpikeProtein = HP_chain(amino_acid_dict,seq[0],True)
	# # print(CoV2SpikeProtein)
    
    # start = time.time()
    # aaSeqCoV = fold(CoV2SpikeProtein)
    # t = (time.time()-start)/60
    # print(f'{len(aaSeqCoV)}aa took {t} min to run')
    # x.append(t)

    # i = [20, 30, 50, 80, 141, 1281]
    i = [20, 30, 50, 80, 141]
    output_graph(i, x, "len(protein)","time(min)","Time vs Protein size-SA", f'time-vs-protein-size-SA({save_num})')

