import os,sys,shutil,glob
import Bio.PDB
import pandas as pd
import subprocess
RES = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
       'THR', 'TRP', 'TYR', 'VAL', 'HIE']


def sheba_run(ref, org):
    cmd3 = 'sheba_01 -x ' + ref + '.pdb ' + org + '.pdb'
    subprocess.call(cmd3, shell=True)
    cmd4 = 'sheba_01 -t ' + org + '.trf ' + org + '.pdb'
    subprocess.call(cmd4, shell=True)
    shutil.copy(org + '.pdb.pdb', org + '_tr.pdb')
    os.remove(org + '.trf')


def ext_chain(aa, bb, cc):
    with open(aa + '_' + bb + '.pdb', 'w') as xx:
        with open(aa + '.pdb', 'r') as xx1:
            lines = xx1.readlines()
            for line in lines:
                if line.startswith('ATOM') > 0 and line[21] == bb:
                    xx.write(line)
        if cc == 'rec':
            xx.write('TER\n')
        elif cc == 'lig':
            xx.write('END\n')


def acc_chem_count(sam, rfeats, folder, ref):
    wdir1 = os.getcwd()

    if not os.path.exists(sam + '_aa.out'):
        os.system('enva.v2 -a ' + sam + '_inp.pdb > ' + sam + '_aa.out')

    mt_raccs = []
    for i in range(len(rfeats)):
        with open(sam + '_aa.out', 'r') as cc:
            lines = cc.readlines()
            for line in lines:
                if line.startswith('ATOM') > 0:
                    aenvs = ' '.join(line[26:].split()).split(' ')
                    if rfeats[i].split('_')[0] == line[22:26].strip() and rfeats[i].split('_')[1] == line[17:20].strip() and \
                            rfeats[i].split('_')[2] == line[12:16].strip():
                        mt_raccs.append(aenvs[5])

    tot_mt_racc = 0
    a_mt_racc = 0
    for mt_racc in mt_raccs:
        a_mt_racc = float(mt_racc) + 1
        tot_mt_racc = tot_mt_racc + a_mt_racc

    os.chdir(folder)

    #	print folder

    idx1 = 0
    pdbs = glob.glob('*_h.pdb')
    atoms = []
    for pdb in pdbs:
        print pdb
        acc = []
        phi = []
        psi = []
        sig = []
        enva_res = ' '
        bbb = pdb.split('.')[0].replace('_h', '').split('_')
        if pdb.find('1ST') > 0 or pdb.find('2ND') > 0:
            enva_res = bbb[0] + 'B' + bbb[1] + '_' + bbb[2] + '.pdb.env'
        else:
            enva_res = bbb[0] + 'B' + bbb[1] + '.pdb.env'
        if not os.path.exists(pdb.split('.')[0].replace('_h', '') + '.out'):
            os.system('enva.v2 -m ' + pdb.split('.')[0].replace('_h', '') + '.pdb > ' + pdb.split('.')[0].replace('_h', '') + '.out')
        if not os.path.exists(pdb.split('.')[0] + '_bb.out'):
            os.system('enva.v2 -b ' + pdb.split('.')[0] + '.pdb > ' + pdb.split('.')[0].replace('_h','') + '_bb.out')
        if not os.path.exists(pdb.split('.')[0].replace('_h', '') + '_aa.out'):
            os.system('enva.v2 -a ' + pdb.split('.')[0].replace('_h', '') + '.pdb > ' + pdb.split('.')[0].replace('_h', '') + '_aa.out')
        if not os.path.exists(pdb.split('.')[0].replace('_h', '') + '.pdb.csv'):
            os.system('python ' + '/lwork02/flexPred/predictFluct.py ' + pdb.split('.')[0].replace('_h', '') + '.pdb XRAY > ' +pdb.split('.')[0].replace('_h', '') + '.pdb.csv')
        if not os.path.exists(enva_res):
            os.system('enva.v2 -e ' + pdb.split('.')[0].replace('_h', '') + '.pdb B')

        raccs = []
        for i in range(len(rfeats)):
            with open(pdb.split('.')[0].replace('_h', '') + '_aa.out', 'r') as af:
                lines = af.readlines()
                for line in lines:
                    if line.startswith('ATOM') > 0:
                        aenvs = ' '.join(line[26:].split()).split(' ')
                        if rfeats[i].split('_')[0] == line[22:26].strip() and rfeats[i].split('_')[1] == line[
                                                                                                         17:20].strip() and \
                                rfeats[i].split('_')[2] == line[12:16].strip():
                            raccs.append(aenvs[5])

        tot_racc = 0
        a_tot_racc = 0
        for racc in raccs:
            a_tot_racc = float(racc) + 1
            tot_racc = tot_racc + a_tot_racc

        del_racc = 0
        del_racc = tot_racc - tot_mt_racc

        tidx = ' '
        tot_lig_acc = 0
        ave_lig_acc = 0
        tot_lig_num = 0
        num = 0
        ratio = 0
        sidx = []
        zdx = 0
        with open(pdb.split('.')[0].replace('_h', '') + '.out', 'r') as ef:
            lines = ef.readlines()
            for line in lines:
                if line.startswith('ATOM') > 0:
                    zdx = 0
                    envs = ' '.join(line[56:].split()).split(' ')
                    if int(envs[8]) == 1:
                        tot_lig_acc = tot_lig_acc + float(envs[0])
                        tot_lig_num = tot_lig_num + 1
                        tidx = envs[3] + '_' + envs[5] + '_' + envs[6]
                        if tidx in rfeats:
                            zdx = rfeats.index(tidx) + 1
                            sidx.append(zdx)
                            num = num + 1

        acc = []
        phi = []
        psi = []
        with open(enva_res, 'r') as ef:
            lines = ef.readlines()
            for line in lines:
                if line.find('chain') < 0 and line.startswith('ATOM') > 0:
                    envs = " ".join(line[56:].split()).split(' ')
                    acc.append(envs[2])
                    phi.append(envs[4])
                    psi.append(envs[5])

        sidx.sort()
        ratio = float(num) / float(len(rfeats))
        if tot_lig_num > 0:
            ave_lig_acc = tot_lig_acc / float(tot_lig_num) + 0.1
        else:
            ave_lig_acc = 0.1

        with open(wdir1 + '/' + sam + '_energy_matrix/total_ac_ct.txt', 'a') as rf1:
            if len(acc) == 9:
                if idx1 == 0:
                    rf1.write('PDB\tP1\tP2\tP3\tP4\tP5\tP6\tP7\tP8\tP9\tPHI1\tPHI2\tPHI3\tPHI4\tPHI5\tPHI6\tPHI7\tPHI8\tPHI9\tPSI1\tPSI2\tPSI3\tPSI4\tPSI5\tPSI6\tPSI7\tPSI8\tPSI9\n')
                rf1.write(
                    '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (pdb.split('.')[0].replace('_h', ''), acc[0], acc[1], acc[2], acc[3], acc[4], acc[5], acc[6], acc[7],acc[8], phi[0], phi[1], phi[2], phi[3], phi[4], phi[5], phi[6], phi[7], phi[8], psi[0], psi[1],psi[2], psi[3], psi[4], psi[5], psi[6], psi[7], psi[8]))
            elif len(acc) == 10:
                if idx1 == 0:
                    rf1.write('PDB\tP1\tP2\tP3\tP4\tP5\tP6\tP7\tP8\tP9\tP10\tPHI1\tPHI2\tPHI3\tPHI4\tPHI5\tPHI6\tPHI7\tPHI8\tPHI9\tPHI10\tPSI1\tPSI2\tPSI3\tPSI4\tPSI5\tPSI6\tPSI7\tPSI8\tPSI9\tPSI10\n')
                rf1.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (pdb.split('.')[0].replace('_h', ''), acc[0], acc[1], acc[2], acc[3], acc[4], acc[5], acc[6], acc[7],acc[8], acc[9], phi[0], phi[1], phi[2], phi[3], phi[4], phi[5], phi[6], phi[7], phi[8], phi[9],psi[0], psi[1], psi[2], psi[3], psi[4], psi[5], psi[6], psi[7], psi[8], psi[9]))

        bfact = []
        with open(pdb.split('.')[0].replace('_h', '') + '.pdb.csv', 'r') as xf:
            lines = xf.readlines()
            for line in lines:
                if line.find('B') > 0:
                    cols = line[:-1].split(',')
                    bfact.append(cols[2])

        with open(wdir1 + '/' + sam + '_energy_matrix/total_bfact.txt', 'a') as bf:
            if len(bfact) == 9:
                if idx1 == 0:
                    bf.write('PDB\tBfactor_P1\tBfactor_P2\tBfactor_P3\tBfactor_P4\tBfactor_P5\tBfactor_P6\tBfactor_P7\tBfactor_P8\tBfactor_P9\n')
                bf.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (pdb.split('.')[0].replace('_h', ''), bfact[0], bfact[1], bfact[2], bfact[3], bfact[4], bfact[5],bfact[6], bfact[7], bfact[8]))
            elif len(bfact) == 10:
                if idx1 == 0:
                    bf.write('PDB\tBfactor_P1\tBfactor_P2\tBfactor_P3\tBfactor_P4\tBfactor_P5\tBfactor_P6\tBfactor_P7\tBfactor_P8\tBfactor_P9\tBfactor_P10\n')
                bf.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (pdb.split('.')[0].replace('_h', ''), bfact[0], bfact[1], bfact[2], bfact[3], bfact[4], bfact[5],bfact[6], bfact[7], bfact[8], bfact[9]))

        with open(wdir1 + '/' + sam + '_energy_matrix/total_rac_ct.txt', 'a') as rf:
            if idx1 == 0:
                rf.write('PDB')
                for rfeat in rfeats:
                    rf.write('\tAA_%s' % (rfeat))
                rf.write('\ttotal_rec_acc\tdelta_racc\tave_lig_acc\t%Match\n')
            for i in range(len(rfeats)):
                if i == 0:
                    rf.write('%s\t%f' % (pdb.split('.')[0].replace('_h', ''), float(raccs[i]) + 1))
                else:
                    rf.write('\t%f' % (float(raccs[i]) + 1))
            rf.write('\t%f\t%f\t%f\t%f\n' % (tot_racc, del_racc, ave_lig_acc, ratio))
        n_bb = 0
        with open(wdir1 + '/' + sam + '_energy_matrix/total_hh_ct.txt', 'a') as hf:
            if idx1 == 0:
                hf.write('PDB\tN.of.BB_full\tN.of.BB_feat\n')
            with open(pdb.split('.')[0].replace('_h', '') + '_bb.out', 'r') as hh:
                hlines = hh.readlines()
                for hline in hlines:
                    for i in range(len(rfeats)):
                        if rfeats[i].split('_')[0] == hline[22:26].strip() and rfeats[i].split('_')[1] == hline[17:20].strip() and \
                                rfeats[i].split('_')[2] == hline[12:16].strip():
                            n_bb = n_bb + 1
            hf.write('%s\t%d\t%d\n' % (pdb.split('.')[0].replace('_h', ''), len(hlines), n_bb))
        idx1 = idx1 + 1
    os.chdir('../')


def conv(aa, cc):
    wdir = os.getcwd()
    res = []
    cres = []
    res1 = ''
    num = []
    cnum = []
    atm = []
    hatm = []
    res2num = {}
    aas = []
    ter = 0
    tpdb = aa.split('_')[0]
    tseq = aa.split('_')[1]
    ttype = aa.split('_')[2]
    aa1 = tpdb + '_' + tseq
    with open(aa1 + '_cl.pdb', 'r') as ff:
        lines = ff.readlines()
        for line in lines:
            if line[12:16].strip() == 'CA' and ter == 0:
                res.append(line[17:20])
                num.append(line[22:26])
            elif line[12:16].strip() == 'CA' and ter > 0:
                aas.append(line[12:16])
                cres.append(line[17:20])
                cnum.append(line[22:26])
            elif line.startswith('TER') > 0:
                ter = ter + 1

    orn = ''
    rn = ''
    idx = 0
    idx1 = 0
    icx = 0
    icx1 = 0
    ter = 0
    cwd = os.getcwd()
    if cc == 'sim_conf':
        conf = aa + '/traj_1/pdb_from_prod'
        os.chdir(conf)
        pdbs = glob.glob('*.pdb.*')
    elif cc == 'PDB_1ST' or cc == 'PDB_2ND':
        conf = aa + '/' + cc
        os.chdir(conf)
        pdbs = glob.glob('*.pdb')
    for pdb in pdbs:
        #	print pdb
        orn = ''
        idx = 0
        idx1 = 0
        icx = 0
        icx1 = 0
        ter = 0
        if cc == 'sim_conf':
            ser = pdb.split('.')[3]
        elif cc == 'PDB_1ST':
            ser = pdb.split('_')[6].split('.')[0] + '_1ST'
        elif cc == 'PDB_2ND':
            ser = pdb.split('_')[6].split('.')[0] + '_2ND'
        with open(aa1 + '-' + ser + '_t.pdb', 'w') as ff1:
            with open(pdb, 'r') as ff2:
                lines = ff2.readlines()
                for line in lines:
                    if line.startswith('ATOM') > 0 and line[17:20] in RES and ter == 0 and line[77] != 'H':
                        if orn == '':
                            TEXT = line[:22] + num[idx] + line[26:]
                            TEXT1 = TEXT[:21] + 'A' + TEXT[22:]
                            ff1.write(TEXT1)
                        elif orn != line[22:26]:
                            idx = idx + 1
                            TEXT = line[:22] + num[idx] + line[26:]
                            TEXT1 = TEXT[:21] + 'A' + TEXT[22:]
                            ff1.write(TEXT1)
                        else:
                            TEXT = line[:22] + num[idx] + line[26:]
                            TEXT1 = TEXT[:21] + 'A' + TEXT[22:]
                            ff1.write(TEXT1)
                        orn = line[22:26]
                    elif line.startswith('ATOM') > 0 and ter > 0 and line[77] != 'H':
                        if orn == '':
                            TEXT = line[:22] + cnum[icx] + line[26:]
                            TEXT1 = TEXT[:21] + 'B' + TEXT[22:]
                            ff1.write(TEXT1)
                        elif orn != line[22:26]:
                            icx = icx + 1
                            TEXT = line[:22] + cnum[icx] + line[26:]
                            TEXT1 = TEXT[:21] + 'B' + TEXT[22:]
                            ff1.write(TEXT1)
                        else:
                            TEXT = line[:22] + cnum[icx] + line[26:]
                            TEXT1 = TEXT[:21] + 'B' + TEXT[22:]
                            ff1.write(TEXT1)
                        orn = line[22:26]
                    elif line.startswith('TER') > 0 or line.startswith('END') > 0:
                        orn = ''
                        icx = 0
                        ter = ter + 1
                        ff1.write(line)

        if cc == 'sim_conf':
            os.system('clean_pdb.py %s-%s_t.pdb A' % (aa1, ser))
            os.system('clean_pdb.py %s-%s_t.pdb B' % (aa1, ser))
            os.system('cat %s-%s_t_A.pdb %s-%s_t_B.pdb > %s-%s.pdb' % (aa1, ser, aa1, ser, aa1, ser))
            os.remove('%s-%s_t_A.pdb' % (aa1, ser))
            os.remove('%s-%s_t_B.pdb' % (aa1, ser))
            os.remove('%s-%s_t_A.fasta' % (aa1, ser))
            os.remove('%s-%s_t_B.fasta' % (aa1, ser))
        else:
            shutil.copy('%s-%s_t.pdb' % (aa1, ser), '%s-%s.pdb' % (aa1, ser))

        with open(aa1 + '-' + ser + '_h.pdb', 'w') as ff3:
            with open(aa1 + '-' + ser + '.pdb', 'r') as ff1:
                lines = ff1.readlines()
                for line in lines:
                    if line.startswith('ATOM') > 0 and line[21] == 'B':
                        TEXT = line.replace('ATOM  ', 'HETATM')
                        ff3.write(TEXT)
                    else:
                        ff3.write(line)
        if cc == 'PDB_1ST' or cc == 'PDB_2ND':
            shutil.copy('%s-%s.pdb' % (aa1, ser), '../dock_res')
            shutil.copy('%s-%s_h.pdb' % (aa1, ser), '../dock_res')
    os.chdir(wdir)


def help():
    print "print help usage\n"
    print "Usage: python %s -s [ sample name ] -i [ input folder] -r [ ref version ] -f [ format ] -d [ cancer or non-cancer]\n" % \
          sys.argv[0]
    return


if len(sys.argv) == 1:
    print "Supply suitable option\n"
    help()
else:

    sam = sys.argv[1]
    # sam = '_'.join([sam1.split('_')[1],sam1.split('_')[0],sam1.split('_')[2]])

    pdb = sam.split('_')[0]
    seq = sam.split('_')[1]
    hlat = sam.split('_')[2]

    temp = pdb + '_' + seq

    try:
        if not os.path.exists(sam + '/' + temp + '_energy_matrix'):
            os.mkdir(sam + '/' + temp + '_energy_matrix', 0777)
    except OSError:
        pass

    try:
        if not os.path.exists(sam + '/dock_res'):
            os.mkdir(sam + '/dock_res')
    except OSError:
        pass

    os.environ['neogear'] = '/lwork02/neoscan_gear/'
    os.environ['rosetta'] = '/lwork02/rosetta_src_2019.40.60963_bundle/'
    #os.environ['flexpred'] = '/lwork02/flexPred/'

    #GEAR = 'neoantigentool/neoscan_gear'
    #ROSETTA_BIN = '/HIP/rosetta_src_2019.40.60963_bundle/tools/protein_tools/scripts'
    #FLEX_PRED = 'neoantigentool/flexPred'
    os.environ['PATH'] += ':' + os.environ['PATH'] + ':' + os.environ['neogear']
    os.environ['PATH'] += ':' + os.environ['PATH'] + ':' + os.environ['rosetta'] + 'main/source/bin/:' + os.environ['rosetta'] + 'tools/protein_tools/scripts/'
    #os.environ['PATH'] += ':' + os.environ['PATH'] + ':' + os.environ['flexpred']
    print(os.environ['PATH'])
    wdir = os.getcwd()

    with open(temp + '_cl.pdb', 'w') as ffx:
        with open(temp + '.pdb', 'r') as ff:
            lines = ff.readlines()
            for line in lines:
                if line.startswith('HETATM') > 0 or line.startswith('ATOM') > 0:
                    if line[16] != ' ':
                        if line[16] == 'A':
                            TEXT = line[:16] + ' ' + line[17:]
                            ffx.write(TEXT)
                    else:
                        ffx.write(line)
                else:
                    ffx.write(line)

    ext_chain(temp + '_cl', 'A', 'rec')
    ext_chain(temp + '_cl', 'B', 'lig')

    if not os.path.exists(pdb + '.out'):
        cmd = 'enva.v2 -m ' + pdb + '.pdb xxxxB > ' + pdb + '.out'
        subprocess.call(cmd, shell=True)

    rfeats = []
    with open(pdb + '.out', 'r') as f1:
        lines = f1.readlines()
        for line in lines:
            if line.startswith('ATOM') > 0 and line[77] != 'H':
                envs = ' '.join(line[56:].split()).split(' ')
                if int(envs[8]) > 0:
                    rfeats.append(envs[3] + '_' + envs[5] + '_' + envs[6])
    print(os.getcwd())
    if os.path.exists(sam + '/traj_1/pdb_from_prod') :
        conv(sam,'sim_conf')
        os.chdir(sam)
        acc_chem_count(sam.replace('_A0201', ''), rfeats, 'traj_1/pdb_from_prod',wdir + '/' + sam.replace('_A0201', '') + '_cl')
        os.chdir(wdir)
        df_rmsd = pd.read_csv(sam + '/' + sam.replace('_A0201', '') + '_energy_matrix/rmsd.txt', sep='\t')
        df_rmsd['PDB'] = df_rmsd['PDB'].str.replace('_B', '')
        df_rac = pd.read_csv(sam + '/' + sam.replace('_A0201', '') + '_energy_matrix/total_rac_ct.txt', sep='\t')
        df_hh = pd.read_csv(sam + '/' + sam.replace('_A0201', '') + '_energy_matrix/total_hh_ct.txt', sep='\t')
        df_ac = pd.read_csv(sam + '/' + sam.replace('_A0201', '') + '_energy_matrix/total_ac_ct.txt', sep='\t')
        df_bfact = pd.read_csv(sam + '/' + sam.replace('_A0201', '') + '_energy_matrix/total_bfact.txt', sep='\t')
        total_df1 = [df_rmsd, df_hh, df_rac, df_ac, df_bfact]
    else:
        conv(sam,'PDB_1ST')
        conv(sam,'PDB_2ND')
        os.chdir(sam)
        acc_chem_count(temp, rfeats, 'dock_res', wdir + '/' + temp + '_cl')
        os.chdir(wdir)
        df_rac = pd.read_csv(sam + '/' + temp + '_energy_matrix/total_rac_ct.txt', sep='\t')
        df_hh = pd.read_csv(sam + '/' + temp + '_energy_matrix/total_hh_ct.txt', sep='\t')
        df_ac = pd.read_csv(sam + '/' + temp + '_energy_matrix/total_ac_ct.txt', sep='\t')
        df_bfact = pd.read_csv(sam + '/' + temp + '_energy_matrix/total_bfact.txt', sep='\t')
        total_df1 = [df_hh, df_rac, df_ac, df_bfact]
    df_final1 = reduce(lambda left, right: pd.merge(left, right, on=['PDB'], how='outer'), total_df1)
    df_final1.to_csv(sam + '/' + temp + '_energy_matrix/' + sam + '_full_env.txt', sep='\t', index=False, na_rep='-')
