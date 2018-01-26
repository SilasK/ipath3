import requests
import matplotlib as mpl
import pandas as pd
import numpy as np
import warnings


values=[-10,-5,0,20]

def map_values_to_colors(values,
                         non_negative_colormap='inferno',
                         divergent_color_map='RdBu_r'):
    """ if values are non negative a sequencial colormap is used, if the values are negative and positive, a symetrical Red to blue colormap is used"""

    if all(np.array(values)>=0):
        norm = None
        cmap= non_negative_colormap
    else:
        extreme= np.abs(values).max()
        norm = mpl.colors.Normalize(vmin=-extreme,vmax=extreme)
        cmap=divergent_color_map

    color_scaler=mpl.cm.ScalarMappable(cmap= cmap,norm=norm)
    rgb=color_scaler.to_rgba(values)[:,:3]
    color_codes= ['RGB({:.0f},{:.0f},{:.0f})'.format(*tuple(row)) for row in rgb*255]

    return color_codes


def map_values_to_range(Values,output_range=[0,1],vmin=None,vmax=None):

    values= np.array(Values)


    if vmin is None: vmin= values.min()
    if vmax is None: vmax= values.max()


    if vmin<0:
        warnings.warn('minimum of values is negative, data will be mapped from [{},{}] -> [{},{}]'.format(vmin,vmax,*tuple(output_range)))

    diff= vmax-vmin



    # clip values
    values[np.where(values>vmax)[0]]=vmax
    values[np.where(values<vmin)[0]]=vmin

    normed_values= (values-vmin)/diff

    assert len(output_range)==2


    out_diff= output_range[1]- output_range[0]

    assert out_diff>0

    return normed_values*out_diff+output_range[0]





def create_selection(data,color_column=None,width_column=None,opacity_column=None,color_kws=None,width_kws=None,opacity_kws=None):
    """transforms a pandas dataframe to selection which can be used in ipath2:


    # Supported data types for indexes
    The following data types are supported by iPath, and can be used to customize the maps. Make sure you use the required prefix for each data type, as shown in the examples below.

    Data type	Prefix	Example ID
    KEGG Pathways	-	00650
    KEGG Compounds	-	C00003
    KEGG KOs	-	K01000
    STRING proteins	-	224324.AQ_626
    KEGG proteins	-	aae:aq_626
    COGs/eggNOGG OGs	-	COG0007
    Enzyme EC numbers	E or EC	E2.4.1.82
    Uniprot IDs/ACCs	UNIPROT:	UNIPROT:Q93015
    IPI IDs	-	IPI00745889
    NCBI GI IDs	ncbi-gi:	ncbi-gi:326314893


    """



    assert type(data)== pd.DataFrame

    output_data= pd.DataFrame(index= data.index)

    if not color_column is None:
        if color_kws is None: color_kws={}
        output_data['color']= map_values_to_colors( data[color_column], **color_kws)

    if not width_column is None:

        width_default_param=dict(output_range=[0,50])
        if not width_kws is None: width_default_param.update(width_kws)

        output_data['Width']= map_values_to_range(data[width_column],**width_default_param).astype(str)
        output_data['Width']="W"+output_data['Width']
    if not opacity_column is None:
        if opacity_kws is None: opacity_kws={}
        output_data['Opacity']= map_values_to_range(data[opacity_column],**opacity_kws).astype(str)

    output_str=''
    for i,row in output_data.iterrows():
        output_str+=' '.join([str(i)]+list(row))+'\n'

    return output_str

## Comunication


def to_parameters(selection ,
                          export_type='svg',
                          include_metabolic=True,
                          include_secondary=False,
                          include_antibiotic=False,
                          include_microbial=False,
                          whole_modules=False,
                          whole_pathways=False,
                          keep_colors=False,
                          default_opacity=1,
                          default_width=3,
                          default_radius=7,
                          default_color='#666666',
                          query_reactions=False,
                          tax_filter='',
                          export_dpi=1200):

    allowed_export_types= ['svg','png','pdf','eps']
    assert export_type in allowed_export_types , "export_type {} needs to be one of {}".froamt(export_type,allowed_export_types)
    #assert map_type=='svg', "I can not save PNG images"

    return dict( selection=selection,
                export_type=export_type,
                keep_colors=int(keep_colors),
                include_metabolic=int(include_metabolic),
                include_secondary= int(include_secondary),
               include_antibiotic= int(include_antibiotic),
               include_microbial= int(include_microbial),
               whole_modules= int(whole_modules),
               default_opacity= default_opacity,
                whole_pathways=int(whole_pathways),
               default_width= default_width,
               default_color= default_color,
               default_radius= default_radius,
               query_reactions= int(query_reactions),
               tax_filter= tax_filter,
               export_dpi=export_dpi)

 # #
def get_map(selection,map_name='map',**param):

    url= 'https://pathways.embl.de/mapping.cgi'



    #print(selection[:300]+'...')
    parameters=to_parameters(selection,**param)

    r = requests.post(url, data=parameters)
    assert r.ok, r.text


    with open(map_name+'.svg','w') as file:
        file.write(r.text)


from IPython.display import SVG

import svgutils.compose as sc

def scale_map(map_name):

    sc.Figure("26cm", "16cm",
        sc.Panel(sc.SVG(map_name+".svg").scale(0.265)),
        ).save(map_name+'_scaled.svg')


def inspect_online(selection,**param):

    url= 'https://pathways.embl.de/ipath.cgi'


    parameters=to_parameters(selection,**param)

    r = requests.post(url, data=parameters)
    return r
