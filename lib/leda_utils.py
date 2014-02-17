from interfits import *


def generateBlankFitsidi(config_xml, n_ant, num_rows, filename_out):
        """ Generate blank FITS-IDI from XML file  """
        
        h1('Creating Primary HDU')
        hdu_primary = make_primary(config=config_xml)
        print hdu_primary.header.ascardlist()
        
        h1('\nCreating ARRAY_GEOMETRY')
        tbl_array_geometry = make_array_geometry(config=config_xml, num_rows=n_ant)
        print tbl_array_geometry.header.ascardlist()

        h1('\nCreating ANTENNA')
        tbl_antenna = make_antenna(config=config_xml, num_rows=self.n_ant)
        print tbl_antenna.header.ascardlist()
        
        h1('\nCreating FREQUENCY')
        tbl_frequency = make_frequency(config=config_xml, num_rows=1)
        print tbl_frequency.header.ascardlist()
        
        h1('\nCreating SOURCE')
        tbl_source = make_source(config=config_xml, num_rows=1)
        print tbl_source.header.ascardlist()
    
        h1('\nCreating UV_DATA')
                # TODO: Fix time and date to julian date
        tbl_uv_data = make_uv_data(config=config_xml, num_rows=num_rows, weights_col=False)
        print tbl_uv_data.header.ascardlist()
        
        h1('Creating HDU list')
        hdulist = pf.HDUList(
                    [hdu_primary, 
                    tbl_array_geometry,
                    tbl_frequency,
                    tbl_antenna,
                    tbl_source, 
                    tbl_uv_data
                    ])
        print hdulist.info()
        
        print('\nVerifying integrity...')            
        hdulist.verify()
  
        if(os.path.isfile(filename_out)):
          print('Removing existing file %s...')%filename_out
          os.remove(filename_out)
        print('Writing to file %s...')%filename_out
        hdulist.writeto(filename_out)

def fillInFits(filename_in, filename_out, uv_pos):
        """ Export data as FITS IDI 
        
        filename_out: str
            output filename
        config_xml: str
            path to config file
        
        """
        
        f_in  = pf.open(filename_in)
        f_out = pf.open(filename_out, mode='update')
        

        tbl_array_geometry_in  = f_in[1]
        tbl_array_geometry_out = f_out[1]
        tbl_antenna_in  = f_in[2]
        tbl_antenna_out = f_out[2]
        tbl_frequency_in  = f_in[3]
        tbl_frequency_out = f_out[3]
        tbl_source_in  = f_in[4]
        tbl_source_out = f_out[4]
        tbl_uv_data_in  = f_in[5]
        tbl_uv_data_out = f_out[5]        
        
        h1('Filling in data')
        h2("ARRAY_GEOMETRY")
        for i in range(len(tbl_array_geometry_in)):
            tbl_array_geometry_out.data[i] = tbl_array_geometry_in.data[i]
        
        h2("ANTENNA")
        for i in range(len(tbl_antenna_in)):
            tbl_antenna_out.data[i] = tbl_antenna_in.data[i]
        
        h2("FREQUENCY")
        for i in range(len(tbl_frequency_in)):
            tbl_frequency_out.data[i] = tbl_frequency_in.data[i]

        h2("SOURCE")
        for i in range(len(tbl_source_in)):
            tbl_source_out.data[i] = tbl_source_out.data[i]
        

        h2("UV_DATA")
        offset = uv_pos * 109
        for i in range(tbl_uv_data_in['DATA'].shape[0]):
            LinePrint("Row %i of %i"%(i+1, self.d_uv_data['DATA'].shape[0]))
            tbl_uv_data_out.data['FLUX'][i] = tbl_uv_data_in.data['FLUX'][i]
            for k in ['UU','VV','WW','BASELINE','DATE']:
                try:
                    tbl_uv_data_out.data[k][i]  = tbl_uv_data_in.data[k][i]
                except:
                    raise
        
        f_out.flush()
        f_out.close()
        f_in.close()