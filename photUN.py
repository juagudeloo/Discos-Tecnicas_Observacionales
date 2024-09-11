import numpy as np

import glob

import os 

from photutils.aperture import CircularAperture
from photutils.aperture import CircularAnnulus
from photutils.aperture import aperture_photometry
from photutils.utils import calc_total_error
from photutils.centroids import centroid_sources, centroid_com

from scipy.stats import mode

from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table
from astropy.stats import sigma_clipped_stats

import pandas as pd

import argparse

import warnings
warnings.filterwarnings("ignore") #for the nan values that are then drop in the final steps.

parser = argparse.ArgumentParser(description='process IRAC (SPITZER) or WISE .fits files to get their photometry.', epilog = 'It is intended to work with MIPS in the future.')


parser.add_argument('instrument', choices = ["IRAC", "WISE"],
                        type=str, help="String type that specifies 'IRAC' or 'WISE' for the instrument to analyze.") 
parser.add_argument('-fd', '--fits_dir', dest='fits_dir', required=True, 
                    type=str, help='A string type that indicates the path to the .fits files.')
parser.add_argument('-cad', '--catalog_dir', dest='catalog_dir', required=True,
                    type=str, help='A string type that indicates the path to the GAIA DR2 catalog.')
parser.add_argument('-od', '--output_dir', dest='output_dir', default="./outputs/",
                    type=str, help='A string type that indicates the path to save the resulting .csv output file.')
parser.add_argument('-r', '--circular_aperture', dest='r', default = 6,
                    type=int, help='Integer value that determines the circular aperture in pixels.')
parser.add_argument('-r_in', '--annular_int_radius', dest='r_in', default = 6,
                    type=int, help='Integer value that determines the annular aperture internal radius in pixels.')
parser.add_argument('-r_out', '--annular_ext_radius', dest='r_out', default = 14,
                    type=int, help='Integer value that determines the annular aperture external radius in pixels.')
parser.add_argument('-v', '--verbose', dest="verbose", action='store_true')

parsed_args = parser.parse_args()

if not parsed_args.fits_dir.endswith("/"):
    parsed_args.fits_dir += "/" 

if not parsed_args.catalog_dir.endswith("/"):
    parsed_args.catalog_dir += "/" 

if not os.path.exists(parsed_args.output_dir):
    os.makedirs(parsed_args.output_dir)

if parsed_args.verbose:
    print(f"""
    Instrument selected:
        > Instrument: {parsed_args.instrument}
        
    Directories specified
        > Fits directory: {parsed_args.fits_dir} 
        > Catalog directory: {parsed_args.catalog_dir} 
        > Output directory: {parsed_args.output_dir} 

    Photometric parameters used
        > Circular aperture radius: {parsed_args.r}
        > Annular aperture internal radius: {parsed_args.r_in}
        > Annular aperture external radius: {parsed_args.r_out}
    """)

def main():    

    #This directory must contain all the files that correspond to one of the instruments (SPITZER, WISE or MIPS).
    fits_dir = parsed_args.fits_dir
    #This directory must contain all the file corresponding to the GAIA catalog (here is used the catalog GAIA DR2 generated using SAO image ds9).
    catalog_dir = parsed_args.catalog_dir

    #########################################
    # Definición de parámetros fotométricos #
    #########################################
    center_box_size = 7
    r = parsed_args.r #Circular aperture in px
    r_in = parsed_args.r_in #Annular aperture external radius
    r_out = parsed_args.r_out #Annular aperture internal radius
    #########################################

    files_list = glob.glob(fits_dir + '*.fits')
    names = []
    for file in files_list:
        if (fits_dir in file) and (file.endswith('.fits')):
            names.append(file.replace(fits_dir,''))

    # Printing the name of the .fits files charged.
    if names != []:
        N_files = len(names)
        print(f'\n Your directory has {N_files} .fits files \n')
        if parsed_args.verbose:
            for l in range(N_files):
                print(f'No. {l+1}: {names[l]}')
    else: 
        print('\n Your directory has not .fits files \n')
                
        print("Names charged in the list!")

    #Charge the file with the catalog generated by SAOimage ds9
    objects = pd.read_csv(catalog_dir+"ds9-ch2.tsv", sep='\t')

    ra = objects["_RAJ2000"]
    dec = objects["_DEJ2000"]
    id_code = objects["Source"]

    decimal_catalog = SkyCoord(ra, dec, unit=(u.degree))
    gaia_catalog = {"_RAJ2000": decimal_catalog.ra.deg,
            "_DEJ2000": decimal_catalog.dec.deg,
            "Source": id_code}
    gaia_catalog = pd.DataFrame(gaia_catalog)

    # Se imprime la tabla en un archivo de texto plano
    all_tables = []
    for k in range(len(files_list)):
        fits_path = files_list[k]
        fits_name = names[k]
        #  catalog = catalog_final

        photom = Photometry_Data_Table(fits_name, fits_path, gaia_catalog, r=r, r_in=r_in, r_out=r_out, center_box_size=center_box_size, verbose=parsed_args.verbose)
        if photom is not None:
            all_tables.append(photom)

    print(f'There are {len(all_tables)} tables from the .fits files')


    #----#  List with all the main objects that are centered by the telescope.
    focus_object = []     
    if parsed_args.instrument == "IRAC":
        identifier = 'OBJECT'
    if parsed_args.instrument == "WISE":    
        identifier = 'COADDID'   

    for m in all_tables:
        ob = m[identifier][0]
        if ob not in focus_object:
            focus_object.append(ob) 
    #----#  Dictionary with each observed object.
    final_filter = {}
    for s in focus_object:
        final_filter[s] = []        # Example: filtro_final = {'SA98':[], 'SA95':[], '[BSA98':[], 'SA101':[], '[ASA98':[], 'SA104':[], 'SA92':[]}

    #----#  Fill the dictionary with the table of each object.
    for n in all_tables:
        for p in focus_object:
            ob = n[identifier][0]
            if ob == p:
                final_filter[ob].append(n.copy())  # Example: filtro_final = {'SA98':[tabla1,tabla2,tabla3,..], ... , 'SA92':[tabla1,tabla2,tabla3,..]}
    
    #----#  Intersection of the objects in the different filters.
    for foc in focus_object:
        current_id = []
        for j in final_filter[foc]:
            current_id.append(j['ID'].data)

        int_d = set(current_id[0]).intersection(*current_id) # Ejemplo para SA98: int_d = {'92_248', ... , '92_347'}

        #----#  Delete the objects that are not in all the filters
        for tab in final_filter[foc]:
            index_of = []
            for i in range(len(tab['ID'])):
                if tab['ID'][i] not in int_d:
                    index_of.append(i)
            tab.remove_rows(index_of)

        #----# Delete the void tables.
        for p in focus_object:
            if len(final_filter[p][0]) == 0:
                del final_filter[p]

        for foc in final_filter.keys():
            let = len(final_filter[foc])
            
            
        #----# Create the tables for each main object.
        for foc in final_filter.keys():
            final_obs_table = pd.DataFrame()
            final_obs_table['OBJECT_ID'] = final_filter[foc][0]['ID']
            final_obs_table['RA'] = final_filter[foc][0]['RA']
            final_obs_table['DEC'] = final_filter[foc][0]['DEC']
            #----# Save the tables as .csv files.
            for j in final_filter[foc]:
                final_obs_table[j.colnames[6]] = j[j.colnames[6]]
                final_obs_table[j.colnames[7]] = j[j.colnames[7]]
                final_obs_table[j.colnames[11] + '_' + j.colnames[6]] = j[j.colnames[11]]
            print(f"Saving .csv for {foc} object")
            final_obs_table.to_csv(parsed_args.output_dir+f'/Table_{foc}_r_{r}_ranul_{r_in}_{r_out}.csv')        


    


def Photometry_Data_Table(fits_name, fits_path, catalog, r, r_in, r_out, center_box_size, verbose, *args):
    """
    Function for
    """
    ################################################
    # Functions
    ################################################
    
    def is_in_pic(w, image, ra, dec):
        """
        Function that selects only the objects from the catalog that are inside the image. 
        """
        ra_max, dec_max = w.array_index_to_world_values(0,0)
        ra_min, dec_min = w.array_index_to_world_values(image.shape[0], image.shape[1])
        if ra_min > ra_max:
            ra_min = w.array_index_to_world_values(0,0)[0]
            ra_max = w.array_index_to_world_values(image.shape[0], image.shape[1])[0]
        if dec_min > dec_max:
            dec_min = w.array_index_to_world_values(0,0)[1]
            dec_max = w.array_index_to_world_values(image.shape[0], image.shape[1])[1]
        return (ra < ra_max) & (ra > ra_min) & (dec < dec_max) & (dec >   dec_min)

    ################################################
    # Main process
    ################################################
    
    with fits.open(fits_path) as FitsData:
        w = WCS(FitsData[0].header)
        image = FitsData[0].data
        fits_header = FitsData[0].header
        if parsed_args.instrument == "IRAC":
            itime  = fits_header['EXPTIME'] 
            ifilter = fits_header['CHNLNUM']  
            DateObs = fits_header['DATE_OBS']
            epadu = fits_header['GAIN']
            #zero points for each channel.
            if ifilter == 1:
                zmag =  18.8027
            elif ifilter == 2:
                zmag = 18.3177
            elif ifilter == 3:
                zmag = 17.8331
            elif ifilter == 4:
                zmag = 17.2120
            else:
                zmag = 17.2120
        if parsed_args.instrument == "WISE":
            ifilter = fits_header['BAND']
            #zero points and electronic gain for each channel.
            if ifilter == 1:
                epadu = 3.20
                itime  = 7.7 #sec
            if ifilter == 2:
                epadu = 3.83
                itime  = 7.7 #sec
            if ifilter == 3:
                epadu = 6.83
                itime = 8.8*fits_header['FRC8P8ET'] \
                + 4.4*fits_header['FRC4P4ET'] \
                + 2.2*fits_header['FRC2P2ET'] \
                + 1.1**fits_header['FRC1P1ET'] #sec
            if ifilter == 4:
                epadu = 24.50
                itime  = 8.8 #sec
            DateObs = fits_header['MIDOBS'] # MID VALUE OF the interval of dates of observation
            #zero points for each channel.
            zmag = fits_header['MAGZP']

    # Writing a txt file that contains the information of the objects that are contained in the picture.
    with open(parsed_args.output_dir+f"Objectlist_{fits_name}.out", "w") as NewListO:
        object_counter = 0 # object counter for the objects that are inside the image.
        for index, row in catalog.iterrows():
            condition = is_in_pic(w, image, row["_RAJ2000"], row["_DEJ2000"])
            if condition:
                object_counter +=1
                X, Y = SkyCoord(row["_RAJ2000"], row["_DEJ2000"], frame="icrs", unit="deg").to_pixel(w)
                NewListO.write(f"{row['_RAJ2000']}    {row['_DEJ2000']}    {row['Source']}   {X}   {Y}   {condition}\n")
        if verbose:
            print(f'\n Found {object_counter} items from the catalog in the file {fits_name} \n')
    if object_counter == 0:
        return None
        quit()

    # Save the coordinates from the objects of the catalog that are in the image.
    with open(parsed_args.output_dir+f"Objectlist_{fits_name}.out", "r") as Obj:
        ListObj = Obj.readlines()

    Final_LO = []
    for obj in ListObj:
        Final_LO.append(obj.split()[:5])
        RA, DEC, ID, x, y = zip(*Final_LO) 
        Final_List = np.array(list(zip(RA,DEC,x,y)), dtype=float)
        ID = np.array(ID,dtype='U20')

    # Drop the objects that are not in the fits files (just in case the is_in_pic function failed numerically)
    mm = [ 0 < i[2] and i[2] < (image.shape[0] - 1) for i in Final_List] # List of [Booleans] (x) where the positions are inside the image.
    ID = ID[mm]
    Final_List = Final_List[mm]
    
    nn = [ 0 < i[3] and i[3] < (image.shape[1] - 1) for i in Final_List] # List of [Booleans] (y) where the positions are inside the image.
    ID = ID[nn]
    Final_List = Final_List[nn]

    # Categorize repeated IDs. 
    u, c = np.unique(ID, return_counts=True)
    dup = u[c > 1]
    for j in dup:
        m = 0
        for i in range(len(ID)):
            if ID[i] == j:
                m += 0.1
                ID[i] = ID[i] + str(m)
    
    # print a preview of the reduced catalog.
    np.set_printoptions(suppress=True, formatter={'float_kind':'{:f}'    .format})

    # Extract the X and Y values.
    _, _, x_init, y_init = zip(*Final_List)

    x, y = centroid_sources(image, x_init, y_init, box_size = center_box_size, centroid_func=centroid_com)
    X, Y = np.array(x), np.array(y)
    NewIDS = np.array(ID) 
        
    # Eliminate the data with NaN or inf centroids.
    is_nan = ~np.isnan(X)
    x, y = X[is_nan],Y[is_nan]      
    Final_List2 = Final_List[is_nan] 
    NewIDS = NewIDS[is_nan]

    # Stars coordinates centroids.
    starloc = list(zip(x,y))
        
    # Each star signal extraction.
    aperture = CircularAperture(starloc, r=r)
    annulus_aperture = CircularAnnulus(starloc, r_in=r_in, r_out=r_out )
    apers = [aperture, annulus_aperture]
    
    # Se genera una tabla de datos.
    phot_table = aperture_photometry(image, apers)

    # The name of the magnitude is assigned depending on the filter from the header.
    name_mag = str(ifilter)

    # Area and flux in the rings. 
    bkg_mean = phot_table['aperture_sum_1'] / annulus_aperture.area
    area_aper = np.array(aperture.area_overlap(image))
    bkg_sum = bkg_mean * area_aper
        
    # Final flux for each object.
    final_sum = phot_table['aperture_sum_0'] - bkg_sum
    phot_table['flux'] = final_sum
    phot_table['flux'].info.format = '%.8g'  
        
    # Instrumental magnitudes.
    phot_table[name_mag + '_mag'] = zmag - 2.5 * np.log10(final_sum) - 2.5 * np.log10(itime)
    phot_table[name_mag + '_mag'].info.format = '%.8g'  
        
    # Instrumental magnitudes error.
    _, _, stdev = sigma_clipped_stats(image, sigma=3.0)
    phot_table[name_mag + '_mag_err'] = 1.0857 * np.sqrt(final_sum/epadu + area_aper*stdev**2 )/final_sum
    phot_table[name_mag + '_mag_err'].info.format = '%.8g'

    # Se agrega a la tabla la RA, DEC, ID y Masa de aire. 
    phot_table['RA'] = [i[0] for i in Final_List2] 
    phot_table['DEC'] = [i[1] for i in Final_List2] 
    phot_table['ID'] = NewIDS
    # phot_table['AIRMASS'] = airmass
    phot_table['DATE-OBS'] = DateObs
    # phot_table['APERTURE'] = r
    # phot_table['Rint'] = r_in
    # phot_table['Rout'] = r_out
    # phot_table['AIRTEMP'] = ccdtemp
    if parsed_args.instrument == "IRAC":
        phot_table['OBJECT'] = fits_header['OBJECT']
    if parsed_args.instrument == "WISE":
        phot_table['COADDID'] = fits_header['COADDID']
    # Se buscan los indices en donde las magnitudes sean NaN y se eliminan
    index_nan = np.argwhere(np.isnan(phot_table[name_mag + '_mag'].data)) 
    phot_table.remove_rows(index_nan)
    filas = len(phot_table)
    return phot_table



if __name__ == "__main__":
    main()