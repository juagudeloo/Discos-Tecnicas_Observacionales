
###################################################################
# Fotometría de imágenes IRAC/SPITZER                             #
# Andrés Felipe Caro Mejía & Giovanni Pinzon                      #
###################################################################

# Librerias
# **photutlis 1.3.0** ; **astropy 4.3.1** 

# INPUT : (1) Directorio en donde se encuentran las imagenes, calibradas
#             astrométricamente lin(27)
#             (https://nova.astrometry.net/upload)
#         (2) Catálogo RA DEC ID "catgaia" lin(51) 
#         (3) Ruta Archivos de Salida *csv conteniendo la fotometría lin (354) 
#         (4) Carga de catalogo landolt lin(359) 
#         (5) Ruta a los archivos *csv lin(375)
#         (6) Parámetros fotométricos lin(254-57)



# OUTPUT : Archivos de texto plano *csv 
           
# INICIO DEL PROGRAMA #
center_box_size=7
           
# Se importa la libreria "glob" con el fin de trabajar con una lista con todos los fits 
import glob
# Directorio en donde se encuentran las imagenes fits  
carpeta = '/media/juanessao2000/hdd/PRINCIPAL-2023-2/MAESTRÍA/Tecnicas_observacionales/Tercer informe - Discos/Data/WISE/'

# Busqueda de los archivos .fits
archivos = glob.glob(carpeta + '*.fits')

# Impresion de los archivos encontados y guardado de nombre del archivo en lista
#print("\n Archivos Encontrados")
nombres = []
for j in archivos:
  if carpeta in j:
    nombres.append(j.replace(carpeta,'')) 

if nombres != []:
  l = len(nombres)
  print(f'\n Su carpeta tiene {l} archivos .fits \n')
  print(f'No. 1: {nombres[0]}')
  print(f'          ....            ')
  print(f'No. {l}: {nombres[l-1]}')
else: 
  print('\n Su carpeta no tiene archivos .fits \n')

# Carga del catálogo con coordenadas e identificador :  RA DEC ID
# !!!!!!! FORMATO DEBE ESTAR EN: RA - DEC - ID   
import numpy as np
archivo_catalogo = "/media/juanessao2000/hdd/PRINCIPAL-2023-2/MAESTRÍA/Tecnicas_observacionales/Tercer informe - Discos/Data/Catalogos_tsv/ds9-W1.tsv" 

catalogo = open(archivo_catalogo,"r")
objects = catalogo.readlines()
catalogo.close()


# Se organiza el catalogo como una lista ordenada
L_O = []
for i in objects:
  L_O.append(i.split())
listObjects = L_O 

# Como el formato del catalogo es hhmmss ggmmss lo pasamos a grados decimales.
ra = [ i[0] for i in L_O ]
dec = [ i[1] for i in L_O ]
id = [ i[2] for i in L_O ]

from astropy.coordinates import SkyCoord
from astropy import units as u
catalogo_decimal = SkyCoord(ra, dec, unit=( u.degree))
catalogo = list(zip(catalogo_decimal.ra.deg,catalogo_decimal.dec.deg,id))
# Fotometría de apertura usando Photutils + Objetos de catálogo
import numpy as np
from photutils.aperture import CircularAperture
from photutils.aperture import CircularAnnulus
from photutils.aperture import aperture_photometry
from photutils.utils import calc_total_error
from scipy.stats import mode
from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table

def Photometry_Data_Table(fits_name, fits_path, catalogo, r, r_in, r_out, center_box_size, *args):
  # Se abre el archivo .fits para guardarlo como una variable i.e. image / fits_data
  F = fits.open(fits_path)
  FitsData = F
  w = WCS(FitsData[0].header)
  fits_data = FitsData[0].data
  fits_header = FitsData[0].header
#print(w.array_shape)
#print(w.world_axis_physical_types)
 # airmass = fits_header['AIRMASS']
  itime  = fits_header['EXPTIME'] 
  ifilter = fits_header['CHNLNUM']  
  DateObs = fits_header['DATE_OBS']
  target = fits_header['OBJECT']
  epadu = fits_header['GAIN']
  F.close()

  image=fits_data

# Funcion que ajusta los objetos del catalogo a solo aquellos que podrían estar en la imágen
  def is_in_pic(w, image, ra, dec):
    ra_max, dec_max = w.array_index_to_world_values(0,0)
    ra_min, dec_min = w.array_index_to_world_values(image.shape[0], image.shape[1])
    if ra_min > ra_max:
      ra_min = w.array_index_to_world_values(0,0)[0]
      ra_max = w.array_index_to_world_values(image.shape[0], image.shape[1])[0]
    if dec_min > dec_max:
      dec_min = w.array_index_to_world_values(0,0)[1]
      dec_max = w.array_index_to_world_values(image.shape[0], image.shape[1])[1]
      
    return (ra < ra_max) & (ra > ra_min) & (dec < dec_max) & (dec >   dec_min)
  NewListO = open(f"Objectlist_{fits_name}.out", "w")
  # Contador de objetos de catálogo que están en la imágen
  object_counter = 0
  for j in range(0,len(catalogo)):
    condicion = is_in_pic(w, image, catalogo[j][0],    catalogo[j][1])
    if condicion:
      object_counter +=1
      X, Y = SkyCoord(catalogo[j][0], catalogo[j][1], frame="icrs", unit="deg").to_pixel(w)
      NewListO.write(f"{catalogo[j][0]}     {catalogo[j][1]}     {catalogo[j][2]}   {X}   {Y}   {condicion}\n")
  NewListO.close()
  print(f'\n Se han encontrado {object_counter} objetos del catalogo en el archivo .fits \n')

  if object_counter == 0:
    return None
    quit()

# Se guardan las coordenadas de los objetos de catálogo que están en la imágen
  Obj = open(f"Objectlist_{fits_name}.out", "r")
  ListObj = Obj.readlines()
  Obj.close()
  Final_LO = []
  for i in ListObj:
    Final_LO.append(i.split()[:5])
  RA, DEC, ID, x, y = zip(*Final_LO) 
  Final_List = np.array(list(zip(RA,DEC,x,y)), dtype=float)
  ID = np.array(ID,dtype='U20')

  # Eliminar los objetos que no esten en el archivo fits (en caso de que la funcion is_in_pic() haya fallado numericamente)
  mm = [ 0 < i[2] and i[2] < (image.shape[0] - 1) for i in Final_List] # Lista de [Booleanos] (x) en las cuales las posiciones si esten en la imagen
  ID = ID[mm]
  Final_List = Final_List[mm]
  nn = [ 0 < i[3] and i[3] < (image.shape[1] - 1) for i in Final_List] # Lista de [Booleanos] (y) en las cuales las posiciones si esten en la imagen
  ID = ID[nn]
  Final_List = Final_List[nn]

# IDs repetidos se categorizan 
  u, c = np.unique(ID, return_counts=True)
  dup = u[c > 1]
  for j in dup:
    m = 0
    for i in range(len(ID)):
      if ID[i] == j:
        m += 0.1
        ID[i] = ID[i] + str(m)

# Se imprime en consola una previsualizacion de como es el nuevo catalogo reducido.
  np.set_printoptions(suppress=True, formatter={'float_kind':'{:f}'    .format})
  print(f"\nSu catalogo reducido es (filas {len(Final_List)}):\n ")
  print("----RA---- ---DEC--- -----x-----  -----y----- -----ID-------\n")
  print(Final_List[0],ID[0])
  print("    ...       ...         ...          ...         ...        \n")
  print(Final_List[len(Final_List)-1],ID[len(Final_List)-1])
  print("---------- --------- ------------  ----------- -------------\n")

# Se extraen los valores X y Y
  _, _, x_init, y_init = zip(*Final_List)




  from photutils.centroids import centroid_sources, centroid_com
  x, y = centroid_sources(image, x_init, y_init, box_size = center_box_size, centroid_func=centroid_com)
  X, Y = np.array(x), np.array(y)
  NewIDS = np.array(ID) 
    
# Se eliminan los datos a los cuales tienen un centroide NaN o inf
  is_nan = ~np.isnan(X)
  x, y = X[is_nan],Y[is_nan]      
  Final_List2 = Final_List[is_nan] 
  NewIDS = NewIDS[is_nan]

# Centroides de coordenadas de las estrellas 
  starloc = list(zip(x,y))
#print(starloc)
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
  print(zmag)  

    
    
# Extracción señal de cada estrella
  aperture = CircularAperture(starloc, r=r)
  annulus_aperture = CircularAnnulus(starloc, r_in=r_in, r_out=r_out )
  apers = [aperture, annulus_aperture]
  # Se genera una tabla de datos.
  phot_table = aperture_photometry(image, apers)

  # Se le asigna nombre de la magnitud dependiendo del filtro en el encabezado
  name_mag = str(ifilter)

  # Area y fujo en los anillos. 
  bkg_mean = phot_table['aperture_sum_1'] / annulus_aperture.area
  area_aper = np.array(aperture.area_overlap(image))
  bkg_sum = bkg_mean * area_aper
    
  # Flujo final para cada objeto
  final_sum = phot_table['aperture_sum_0'] - bkg_sum
  phot_table['flux'] = final_sum
  phot_table['flux'].info.format = '%.8g'  
    
  # Magnitudes Instrumentales
  phot_table[name_mag + '_mag'] = zmag - 2.5 * np.log10(final_sum) + 2.5 * np.log10(itime)
  phot_table[name_mag + '_mag'].info.format = '%.8g'  
    
  # Error de las Magnitudes Instrumentales
  from astropy.stats import sigma_clipped_stats
  mean, median, std = sigma_clipped_stats(image, sigma=3.0)
  stdev = std
  phot_table[name_mag + '_mag_err'] = 1.0857 * np.sqrt(final_sum    /epadu + area_aper*stdev**2 )/final_sum
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
  phot_table['OBJECT'] = fits_header['OBJECT']
  # Se buscan los indices en donde las magnitudes sean NaN y se eliminan
  index_nan = np.argwhere(np.isnan(phot_table[name_mag + '_mag'].data)) 
  phot_table.remove_rows(index_nan)
  filas = len(phot_table)
  return phot_table


#########################################
# Definición de parámetros fotométricos #
#########################################
r = 6 #Apertura en px
r_in = 6 #Radio interno anillo cielo
r_out = 14 #Radio externo
#########################################
# Se imprime la tabla en un archivo de texto plano
all_tables = []
for k in range(len(archivos)):
  fits_path = archivos[k]
  fits_name = nombres[k]
#  catalogo = catalogo_final

  photom = Photometry_Data_Table(fits_name, fits_path, catalogo, r=r, r_in=r_in, r_out=r_out, center_box_size=center_box_size)
  if photom is not None:
    all_tables.append(photom)

print(f'Se tienen {len(all_tables)} tablas de las imagenes .fits')
#----# Crea lista con los nombres de los objetos a los cuales se enfoca el telescopio
focus_object = []             
for m in all_tables:
  
 #if m != []:
  ob = m['OBJECT'][0]
  
  if ob not in focus_object:
    
    focus_object.append(ob)   # Ejemplo: focus_object = ['SA98', 'SA95', '[BSA98', 'SA101', '[ASA98', 'SA104', 'SA92']

#----#  Se crea diccionario con cada objeto de enfoque
filtro_final = {}
for s in focus_object:
  filtro_final[s] = []        # Ejemplo: filtro_final = {'SA98':[], 'SA95':[], '[BSA98':[], 'SA101':[], '[ASA98':[], 'SA104':[], 'SA92':[]}

#----#  Se llena el diccionario
for n in all_tables:
  for p in focus_object:
    ob = n['OBJECT'][0]
    if ob == p:
      filtro_final[ob].append(n.copy())  # Ejemplo: filtro_final = {'SA98':[tabla1,tabla2,tabla3,..], ... , 'SA92':[tabla1,tabla2,tabla3,..]}

#----#  Para cada observacion de enfoque se hace la interseccion de los objetos que esten en los tres filtros
for foc in focus_object:
  current_id = []
  for j in filtro_final[foc]:
    current_id.append(j['ID'].data)
  
  int_d = set(current_id[0]).intersection(*current_id) # Ejemplo para SA98: int_d = {'92_248', ... , '92_347'}

  #----# Se eliminan los objetos que no esten en los tres filtros
  for tab in filtro_final[foc]:
    index_of = []
    for i in range(len(tab['ID'])):
      if tab['ID'][i] not in int_d:
        index_of.append(i)
    tab.remove_rows(index_of)

#----# Eliminar las tablas que esten vacias
for p in focus_object:
  if len(filtro_final[p][0]) == 0:
    del filtro_final[p]

for foc in filtro_final.keys():
  let = len(filtro_final[foc])
  
    
from astropy.table.table import QTable

# noche = input('Ingrese la noche en que hizo la fotometria')
noche = 'N3'
''
#----# Se crean tablas para cada objeto de enfoque
for foc in filtro_final.keys():
  final_obs_table = QTable()
  final_obs_table['OBJECT_ID'] = filtro_final[foc][0]['ID']
  final_obs_table['RA'] = filtro_final[foc][0]['RA']
  final_obs_table['DEC'] = filtro_final[foc][0]['DEC']

#----# Se guardan las tablas como archivos .csv
  counter = 0
  for j in filtro_final[foc]:
    
  
    final_obs_table[j.colnames[6] + '_' + str(counter//3)] = j[j.colnames[6]]
    final_obs_table[j.colnames[7] + '_' + str(counter//3)] = j[j.colnames[7]]
    final_obs_table[j.colnames[11] + '_' + j.colnames[6] + '_' + str(counter//3)] = j[j.colnames[11]]
    counter += 1
  final_obs_table.write(f'/media/juanessao2000/hdd/PRINCIPAL-2023-2/MAESTRÍA/Tecnicas_observacionales/Tercer informe - Discos/Data/WISE_outputs/Table_{foc}.csv', overwrite=True)    




# FIN DEL PROGRAMA #
# Observatorio Astronómico Nacional 2024 #
