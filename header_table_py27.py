from astropy.table import Table, Column
from astropy.io import fits
import argparse
import os
import numpy as np
import glob
from multiprocessing import Pool


def Staistics_table_main(files):
	file_name     = []
	file_date     = []
	file_aperture = []
	file_postarg1 = []
	file_postarg2 = []
	file_pa_v3    = []
	file_proposal = []
	file_ra_targ  = []
	file_dec_targ = []
	file_bunit    = []
	file_orientat = []
	file_exptm    = []
	file_filter   = []



	for im in files:
		hdu = fits.open(im)
		date = hdu[0].header['date-obs']
		filter_name=hdu[0].header['filter']
		rootname = hdu[0].header['rootname']
		exptm    = hdu[0].header['exptime']
		postarg1 = hdu[0].header['POSTARG1']
		postarg2 = hdu[0].header['POSTARG2']
		pa_v3    = hdu[0].header['PA_V3']
		aperture = hdu[0].header['APERTURE']
		proposal = hdu[0].header['PROPOSID']
		ra_targ  = hdu[0].header['RA_TARG']
		dec_targ = hdu[0].header['DEC_TARG']
		bunit    = hdu[1].header['BUNIT']
		orientat = hdu[1].header['ORIENTAT']

		file_name     = np.append(file_name,rootname)
		file_filter   = np.append(file_filter,filter_name)
		file_date     = np.append(file_date,date)
		file_exptm    = np.append(file_exptm,exptm)
		file_aperture = np.append(file_aperture,aperture)
		file_postarg1 = np.append(file_postarg1, postarg1) 
		file_postarg2 = np.append(file_postarg2, postarg2) 
		file_pa_v3    = np.append(file_pa_v3   , pa_v3) 
		file_proposal = np.append(file_proposal, proposal) 
		file_ra_targ  = np.append(file_ra_targ , ra_targ) 
		file_dec_targ = np.append(file_dec_targ, dec_targ) 
		file_bunit    = np.append(file_bunit   , bunit) 
		file_orientat = np.append(file_orientat, orientat)  

		



	
	t1=Table()
	t1['Name'] = file_name
	t1['Date'] = file_date
	t1['ProposalID'] = file_proposal
	t1['Exposure Time'] = file_exptm
	t1['Filter'] = file_filter
	t1['RA'] = file_ra_targ
	t1['Dec'] = file_dec_targ
	t1['PA_V3'] = file_pa_v3
	t1['Postarg1'] = file_postarg1
	t1['Postarg2'] = file_postarg2
	t1['Orientat'] = file_orientat
	t1['Aperture'] = file_aperture
	t1.write('{}_ext4_flt_statistic_table.txt'.format(filter_name),format='ascii.fixed_width')

#	t2=Table()
#	t2['Name'] = file_name
#	t2['Date'] = file_date
#	t2['Filter'] = file_filter
#	t2['Aperture'] = file_aperture
#	t2['Exposure Time'] = file_exptm
#	t2['Max'] = file_max2
#	t2['Min'] = file_min2
#	t2['Mean'] = file_mean2
#	t2['STDDEV'] = file_STDDEV2
#	t2.write('{}_ext1_flt_statistic_table.txt'.format(filter_name1),format='ascii.fixed_width', overwrite=True)




def parse_args():
    """Parses command line arguments.

    Parameters:
        nothing

    Returns:
        args : argparse.Namespace object
            An argparse object containing all of the added arguments.

    Outputs:
        nothing
    """

    #Create help string:
    path_help = 'Path to the folder with files to run tweakreg.'
    # Add arguments:
    parser = argparse.ArgumentParser()
    parser.add_argument('--path', '-path', dest = 'path', action = 'store',
                        type = str, required = True, help = path_help)


    # Parse args:
    args = parser.parse_args()

    return args

if __name__ == '__main__':
    args = parse_args()
    path=args.path
    os.chdir(path)
    files = sorted(glob.glob('*.fits'))
    Staistics_table_main(files)
    os.system('mv *.txt ../')
    os.system('cd ../')

