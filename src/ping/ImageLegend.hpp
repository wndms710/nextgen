#ifndef IMAGELEGEND_HPP
#define IMAGELEGEND_HPP


#include <stdio.h>
#include <string.h>
#include <stdlib.h>


#ifdef WITH_FREETYPE
  #include <ft2build.h>
  #include FT_FREETYPE_H
  #include FT_GLYPH_H
  static const string legend_font_file = WITH_FREETYPE;
#else
  #include "CtiType.hpp"
  static const string legend_font_file = "";
#endif

static bool b_legend_font = true;

class ImageLegend {

private:

  bool b_init;


  int image_nx; //full image dimensions
  int image_ny; //full image dimensions

  int lg_off_x; //offset in pixels of left of legend
  int lg_off_y; //offset in pixels of top of legend

  float lg_ar; //lg_nx/lg_ny
  float lg_sc; //lg_ny/ny

  int lg_min_y;  //min height of legend in pixels
  int lg_min_x;  //min width of legend in pixels

  int lg_ny;  //full height in pixels of legend
  int lg_nx;  //full width in pixels of legend

  int lg_txt; //number of pixels for text above, below, and left of colorbar
  int lg_bar; //pixel shift of colorbar

  //buffer, overlay on image using transparency
  unsigned char (*rgba)[4];

public:

  ImageLegend() {
    b_init = false;

    lg_off_x = 10;
    lg_off_y = 10;
    
    lg_ar = 0.50; //lg_nx/lg_ny
    lg_sc = 0.25; //lg_ny/ny
    
    lg_min_y = 200;
    lg_min_x = lg_min_y*lg_ar;

    rgba = NULL;
  }

  ~ImageLegend() {
    if (rgba){
      delete[] rgba;
      rgba = NULL;
    }
  }

  void init(const int _image_nx, const int _image_ny, const double offset_frac_x, const double offset_frac_y, const ColorMap &theColorMap){
    init(_image_nx, _image_ny, theColorMap);
//    int off_y_max = (int) ((1.0-lg_sc) * image_ny - 10);
//    int off_x_max = (int) (image_nx - (lg_ar * (lg_sc * image_ny)) - 10);
    int off_y_max = (int) (image_ny - lg_ny - 10);
    int off_x_max = (int) (image_nx - lg_nx - 10);

    lg_off_x = std::min(std::max(10,(int) (image_nx * offset_frac_x)),off_x_max);
    lg_off_y = std::min(std::max(10,(int) (image_ny * offset_frac_y)),off_y_max);
  }

  //defaults
  void init(const int _image_nx, const int _image_ny, const ColorMap &theColorMap){ 
    b_init = true;

    image_ny = _image_ny;
    image_nx = _image_nx;
 
    //TODO: warn here image is too small for a legend
    assert(image_nx > lg_min_x && image_ny > lg_min_y);

    lg_ny = (int) (lg_sc * image_ny);
    if (lg_ny < lg_min_y)
      lg_ny = lg_min_y;

    lg_nx = (int) (lg_ar * lg_ny);   
    lg_txt = (int) (0.09 * lg_ny);
    lg_bar = (int) (0.20 * lg_nx);

    rgba = new unsigned char[lg_nx*lg_ny][4];
    for (int ipx=0; ipx<lg_nx*lg_ny; ++ipx){
      rgba[ipx][0] = 73;
      rgba[ipx][1] = 175;
      rgba[ipx][2] = 205;
      rgba[ipx][3] = 0;
    }

    float cbar_height_inv = 1.0f/(lg_ny - lg_txt - lg_txt);

    for (int irow=lg_txt; irow<(lg_ny-lg_txt); ++irow){
      float phi = 1.0f - (irow-lg_txt)*cbar_height_inv; //0 to 1
      unsigned char my_rgb[3];
      theColorMap.calcColorMapRGB(my_rgb, phi);
      for (int icol=lg_bar; icol<(lg_nx-lg_bar); ++icol){
         int ipx = irow*lg_nx + icol;
         FOR_I3 rgba[ipx][i] = my_rgb[i];
         rgba[ipx][3] = 255;
      }
    }

  }

  bool isInit() {
   return b_init;
  }

  void addLegend(unsigned char (*image_rgb)[3], PngDataChunk * zoNe, const float rmin, const float rmax, const string &varName){
//    cout << "Adding Legend..." << endl;
  
    ostringstream ss;
    if (rmin < 0.0f || rmax< 0.0f) ss << std::showpos;

    ss << std::scientific << std::setprecision(3) << rmin;
    string conf_text_min = ss.str();

    ss.str(""); 
//ss.clear();
//    ss << std::scientific << std::setprecision(3) << rmax;
    ss << rmax;
    string conf_text_max = ss.str();

    
    if (b_legend_font) renderString(conf_text_max, 0, 0);
    if (b_legend_font) renderString(conf_text_min, 0, lg_ny-lg_txt);
    if (b_legend_font){
      if (varName.length()>3)
        if (varName.length()>10)
          renderString(varName.substr(0,10),(int)(0.5*(lg_nx-lg_txt)),lg_ny-lg_txt,true);
        else
          renderString(varName,(int)(0.5*(lg_nx-lg_txt)),lg_ny-lg_txt,true);
      else
        renderString(varName,0,(int) (0.5*(lg_ny-lg_txt)));
    }


    
    unsigned char zoNe_uc[2] = {(unsigned char) 255, (unsigned char) 255};
    for (int irow = 0; irow<lg_ny; ++irow){
      int irow_image = irow + lg_off_y;
      for (int icol = 0; icol<lg_nx; ++icol){
         int icol_image = icol + lg_off_x;
         int ipx = irow*lg_nx + icol;
         int ipx_image = irow_image*image_nx + icol_image;
         float alpha = rgba[ipx][3]/255.0f;
         if (alpha>0){
           FOR_I3 image_rgb[ipx_image][i] = std::max(std::min((int) (alpha * rgba[ipx][i] + (1.0f-alpha)*image_rgb[ipx_image][i]),255),0);
           zoNe->set(ipx_image,zoNe_uc);
         }
      }
    }
  }

private:


#ifndef WITH_FREETYPE
//  int renderString(const string &theText, const int i0, const int j0, const bool b_rotateAndCenter=false){
//    cout << "*******************************************************************************" << endl;
//    cout << "Warning: no text will be displayed on image legend" << endl;
//    cout << "  To write a legend with text (i.e. variable ranges) compile with libfreetype" << endl;
//    cout << "  For more information see https://www.freetype.org/" << endl;
//    cout << "  Set a flag and font file at compile time in Makefile.in:" << endl;
//    cout << "    -DWITH_FREETYPE='\"/usr/share/fonts/truetype/freefont/FreeSans.ttf\"'" << endl;
//    cout << "*******************************************************************************" << endl;
//    b_legend_font = false;
//    return 0;
//  }
//
    int renderString(const string &theText, const int i0, const int j0, const bool b_rotateAndCenter=false){
      //cti_type_height = 18  200*0.09 = 18

      if (cti_type_height[0] > lg_txt){
        cout << "*******************************************************************************" << endl;
        cout << "Warning: no text will be displayed on image legend." << endl;
        cout << "  Character map to large for legend box." << endl;
        cout << "*******************************************************************************" << endl;
        b_legend_font = false;
        return 0;
      }

      //lg_txt is the height we have to work with, look in our font array for the corresponding char_width
      int ifs;
      int iscale = 0;
      int sizeFlag = -1;
      while (sizeFlag<0){
        ++iscale;
        for (ifs=0; ifs<cti_type_n; ++ifs){
          if (cti_type_height[ifs]*iscale==lg_txt){
            sizeFlag=0;
            break;
          }
          else if (cti_type_height[ifs]*iscale>lg_txt){
            sizeFlag=1;
            break;
          }
        }
      }
      if (sizeFlag==1){
        //cout << "*******************************************************************************" << endl;
        //cout << "Warning: requested font size " << lg_txt << "px using " << iscale*cti_type_height[ifs-1] << "px" << endl;
        //cout << "*******************************************************************************" << endl;
        --ifs;
      }
      //if (iscale>1)
      //  cout << "Scaling font size " <<cti_type_height[ifs]<< "px to " << cti_type_height[ifs]*iscale<< "px" << endl;

      int char_width = cti_type_width[ifs]*iscale;
      int char_height= cti_type_height[ifs]*iscale;
      //double scale = ((double) lg_txt)/cti_type_height;
      //int char_width = scale*cti_type_width;

      //cout << " scale " << scale << " char_width " << char_width << " char_height " << lg_txt << endl;

      int i,j;
      if (b_rotateAndCenter){
        i = i0;
        j = j0 - 0.5*(lg_ny-2.0*lg_txt) + 0.5*char_width*theText.length();
      }
      else{
        i = i0 + 0.5*lg_nx - 0.5*char_width*theText.length();
        j = j0;
      }
      for (int ii = 0; ii < theText.length(); ++ii) {
        const int ichar = int(theText[ii])-32;
        assert((ichar >= 0)&&(ichar < cti_type_count));
        for (int di = 0; di < char_width; ++di) {
          for (int dj = 0; dj < char_height; ++dj) {
            int ipx;
            if (b_rotateAndCenter)
              ipx = (j-di)*lg_nx + (i+dj);
            else
              ipx = (j+dj)*lg_nx + (i+di);

            //index into cti_type_data; 
            int di_type = di/iscale;  
            int dj_type = dj/iscale; 

            float alpha0 = (255 - (int) cti_type_data_ptr[ifs][ichar*cti_type_width[ifs]*cti_type_height[ifs]+dj_type*cti_type_width[ifs]+di_type])/255.0f;
            float alpha1 = rgba[ipx][3]/255.0f;
            FOR_I3 rgba[ipx][i] = std::max(std::min((int) ( ((1.0f-alpha0)*alpha1*rgba[ipx][i])/(alpha0 + alpha1*(1.0f-alpha0))  ),255),0);
            rgba[ipx][3] =  std::max(std::min((int) ( (alpha0 + alpha1*(1.0f-alpha0))*255.0f ),255),0);
          }
        }
        if (b_rotateAndCenter)
          j -= char_width;
        else
          i += char_width;
      }
      return 1;
    }

#else

  int renderString(const string &theText, const int i0, const int j0, const bool b_rotateAndCenter=false){
    const int MAX_GLYPHS = 20; //TODO
     
    const int conf_size = lg_txt;
    //string conf_font_file = "/usr/share/fonts/truetype/freefont/FreeSerif.ttf";
    //string conf_font_file = "/usr/share/fonts/truetype/freefont/FreeMono.ttf";
    //string conf_font_file = "/usr/share/fonts/truetype/freefont/FreeSans.ttf";


    const int num_chars =  theText.length();

    FT_Face face; 
    FT_Library ft;
    FT_Error err;

    FT_Init_FreeType(&ft);
    err = FT_New_Face(ft, legend_font_file.c_str(), 0, &face);
    if (err) {
      cout << "*******************************************************************************" << endl;
      cout << "Warning: unable to load font file for legend, no text will be displayed" << endl;
      cout << "         Set font file at compile time in Makefile.in i.e." << endl;
      cout << "         -DWITH_FREETYPE='\"/usr/share/fonts/truetype/freefont/FreeSans.ttf\"'" << endl;
      cout << "*******************************************************************************" << endl;
      b_legend_font = false;
      return 0;
    }
    FT_Set_Pixel_Sizes(face, 0, conf_size);
 
    FT_GlyphSlot  slot = face->glyph;   /* a small shortcut */
    FT_UInt       glyph_index;
    FT_Bool       use_kerning;
    FT_UInt       previous;
    int           pen_x, pen_y, n;
    
    FT_Glyph      glyphs[MAX_GLYPHS];   /* glyph image    */
    FT_Vector     pos   [MAX_GLYPHS];   /* glyph position */
    FT_UInt       num_glyphs;
  
    pen_x = 0;   /* start at (0,0) */
    pen_y = 0;
    
    num_glyphs  = 0;
    use_kerning = FT_HAS_KERNING( face );
    previous    = 0;
    
    for ( n = 0; n < num_chars; n++ ) {
      /* convert character code to glyph index */
      glyph_index = FT_Get_Char_Index( face, theText[n] );
    
      /* retrieve kerning distance and move pen position */
      if ( use_kerning && previous && glyph_index )
      {
        FT_Vector  delta;
    
    
        FT_Get_Kerning( face, previous, glyph_index,
                        FT_KERNING_DEFAULT, &delta );
    
        pen_x += delta.x >> 6;
      }
    
      /* store current pen position */
      pos[num_glyphs].x = pen_x;
      pos[num_glyphs].y = pen_y;
    
      /* load glyph image into the slot without rendering */
      err = FT_Load_Glyph( face, glyph_index, FT_LOAD_DEFAULT );
      if ( err )
        continue;  /* ignore errors, jump to next glyph */
    
      /* extract glyph image and store it in our table */
      err = FT_Get_Glyph( face->glyph, &glyphs[num_glyphs] );
      if ( err )
        continue;  /* ignore errors, jump to next glyph */
    
      /* increment pen position */
      pen_x += slot->advance.x >> 6;
    
      /* record current glyph index */
      previous = glyph_index;
    
      /* increment number of glyphs */
      num_glyphs++;
    }

    /* compute string dimensions in integer pixels */
    FT_BBox string_bbox;
    compute_string_bbox( &string_bbox, pos, glyphs, num_glyphs );

    int string_width  = string_bbox.xMax - string_bbox.xMin;
    int string_height = string_bbox.yMax - string_bbox.yMin;
  
    /* compute start pen position in pixels */
     
    int start_x, start_y;
  
    if (!b_rotateAndCenter){
      start_x = (int) (( 0.5*( lg_nx - string_width ) ) * 64);
      start_y = (int) (( 0.5*( lg_txt - string_height ) ) * 64);
    }
    else{
      start_x = (int) (( 0.5*( lg_ny - 2*lg_txt  - string_width  ) ) * 64);
      start_y = (int) (( 0.5*( lg_txt - string_height ) ) * 64); //DAPDAP
    }

    for ( n = 0; n < num_glyphs; n++ )
    {
      FT_Glyph   image;
      FT_Vector  pen;
    
      image = glyphs[n];
    
      pen.x = start_x + (pos[n].x << 6);
      pen.y = start_y + (pos[n].y << 6);

//      cout << "Pen x y: " << pen.x << " " << pen.y << endl; //DAP
    
      err = FT_Glyph_To_Bitmap( &image, FT_RENDER_MODE_NORMAL,
                                &pen, 0 );
      if ( !err )
      {
        FT_BitmapGlyph  bit = (FT_BitmapGlyph)image;

//        cout << "bit->left " << bit->left << " bit->top " << bit->top << " bit->bitmap.rows " << bit->bitmap.rows << " bit->bitmap.width " << bit->bitmap.width << " bit->bitmap.pitch " << bit->bitmap.pitch << endl; //DAP
    
        for (int irow = 0; irow < bit->bitmap.rows; ++irow) {
          for (int icol = 0; icol < bit->bitmap.width; ++icol) {
            const unsigned char pixel = bit->bitmap.buffer[irow*bit->bitmap.pitch + icol];
              int rgba_row, rgba_col;
           
            if (b_rotateAndCenter){
              rgba_row =  j0 - (bit->left + icol);
              rgba_col =  i0 + (lg_txt - bit->top) + irow;
            }
            else{
              rgba_col =  i0 + bit->left + icol;
              rgba_row =  j0 + (lg_txt - bit->top) + irow;
            }

            int ipx = rgba_row*lg_nx + rgba_col;
            if (rgba_col<lg_nx){
              //FOR_I3 rgba[ipx][i] = 0; //black
              //rgba[ipx][3] = pixel;
              float alpha0 = pixel/255.0f;
              float alpha1 = rgba[ipx][3]/255.0f;
              FOR_I3 rgba[ipx][i] = std::max(std::min((int) ( ((1.0f-alpha0)*alpha1*rgba[ipx][i])/(alpha0 + alpha1*(1.0f-alpha0))  ),255),0);
              rgba[ipx][3] =  std::max(std::min((int) ( (alpha0 + alpha1*(1.0f-alpha0))*255.0f ),255),0);
            }
            else
              cout << "Warning: font pixel out of range" << endl; //TODO
          }
        }
    
        FT_Done_Glyph( image );
        FT_Done_Glyph( glyphs[n] );
      }
    }

    FT_Done_FreeType(ft);
    return 1;
  }

  const char *render_glyph(FT_Face *face, const int conf_size, const char * conf_font_file, const char * theText, bool conf_anti_alias) {

    FT_Library ft;
    FT_Error err;
   
    err = FT_Init_FreeType(&ft);
    if (err) return "freetype init error";
   
    err = FT_New_Face(ft, conf_font_file, 0, face);
    if (err == FT_Err_Unknown_File_Format)
      return "unknown font file format";
    else if (err)
      return "error reading font file";
   
    err = FT_Set_Pixel_Sizes(*face, 0, conf_size);
    if (err) return "error setting font size";
   
    FT_UInt index = FT_Get_Char_Index(*face, *theText);
    if (index == 0) return "no glyph found for char";
   
    err = FT_Load_Glyph(*face, index, FT_LOAD_DEFAULT);
    if (err) return "error loading glyph";
   
    err = FT_Render_Glyph((*face)->glyph, conf_anti_alias ?
                          FT_RENDER_MODE_NORMAL :
                          FT_RENDER_MODE_MONO);
    if (err) return "error rendering glyph";
   
    return NULL;
  }
 

  void  compute_string_bbox( FT_BBox  *abbox, FT_Vector *pos, FT_Glyph *glyphs, FT_UInt num_glyphs ) {
    FT_BBox  bbox;
    FT_BBox  glyph_bbox;
  
    /* initialize string bbox to "empty" values */
    bbox.xMin = bbox.yMin =  32000;
    bbox.xMax = bbox.yMax = -32000;
  
    /* for each glyph image, compute its bounding box, */
    /* translate it, and grow the string bbox          */
    for (int n = 0; n < num_glyphs; n++ )
    {
      FT_Glyph_Get_CBox( glyphs[n], ft_glyph_bbox_pixels,
                         &glyph_bbox );
  
      glyph_bbox.xMin += pos[n].x;
      glyph_bbox.xMax += pos[n].x;
      glyph_bbox.yMin += pos[n].y;
      glyph_bbox.yMax += pos[n].y;
  
      if ( glyph_bbox.xMin < bbox.xMin )
        bbox.xMin = glyph_bbox.xMin;
  
      if ( glyph_bbox.yMin < bbox.yMin )
        bbox.yMin = glyph_bbox.yMin;
  
      if ( glyph_bbox.xMax > bbox.xMax )
        bbox.xMax = glyph_bbox.xMax;
  
      if ( glyph_bbox.yMax > bbox.yMax )
        bbox.yMax = glyph_bbox.yMax;

    }
  
    /* check that we really grew the string bbox */
    if ( bbox.xMin > bbox.xMax )
    {
      bbox.xMin = 0;
      bbox.yMin = 0;
      bbox.xMax = 0;
      bbox.yMax = 0;
    }
  
    /* return string bbox */
    *abbox = bbox;

  }
#endif


};

#endif

