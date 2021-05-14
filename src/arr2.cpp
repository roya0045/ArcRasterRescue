#include <zlib.h>

#include <bitset>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>


#pragma once

#include <gdal_priv.h>
#include "cpl_string.h"
#include "gdal_frmts.h"
#include <algorithm>
#include <cassert>
#include <chrono>
#include <fstream>
#include <typeinfo>



struct RasterFields {
  double raster_mtolerance;
  double raster_xytolerance;
  double raster_zorig;
  double raster_morig;
  double raster_mscale;
  double raster_zscale;
  double raster_xorig;
  double raster_yorig;
  double raster_xyscale;
  double raster_ztolerance;
  bool   raster_has_m;
  bool   raster_has_z;
  std::string wkt;
  std::string raster_column;
};

struct Shape {
  double ymax;
  double xmax;
  double xmin;
  double ymin;
  double morig;
  double zorig;
  double zscale;
  double mscale;
  double xyscale;
  double xorig;
  double yorig;
  bool has_z;
  bool has_m;
  double mtolerance;
  double ztolerance;
  double xytolerance;
  std::string wkt;
};

struct Field {
  std::string  name;
  std::string  alias;
  int8_t       type;
  bool         nullable;
  RasterFields raster;
  Shape        shape;
  void print() const;
};


/* ==================================================================== */
/************************************************************************/

class JPEG2000Dataset final: public GDALJP2AbstractDataset
{
    friend class JPEG2000RasterBand;

    jas_stream_t *psStream;
    jas_image_t *psImage;
    int         iFormat;
    int         bPromoteTo8Bit;

    int         bAlreadyDecoded;
    int         DecodeImage();

  public:
                JPEG2000Dataset();
                ~JPEG2000Dataset();

    static int           Identify( GDALOpenInfo * );
    static GDALDataset  *Open( GDALOpenInfo * );
};

/************************************************************************/
/* ==================================================================== */
/*                            JPEG2000RasterBand                        */
/* ==================================================================== */
/************************************************************************/

class JPEG2000RasterBand final: public GDALPamRasterBand
{
    friend class JPEG2000Dataset;

    // NOTE: poDS may be altered for NITF/JPEG2000 files!
    JPEG2000Dataset     *poGDS;

    jas_matrix_t        *psMatrix;

    int                  iDepth;
    int                  bSignedness;

  public:

    JPEG2000RasterBand( JPEG2000Dataset *, int, int, int );
    virtual ~JPEG2000RasterBand();

    virtual CPLErr IReadBlock( int, int, void * ) override;
    virtual GDALColorInterp GetColorInterpretation() override;
};



class GDBTable {
 public:
  std::ifstream gdbtable, gdbtablx;
  int32_t nfeaturesx;
  int32_t size_tablx_offsets;
  std::vector<Field> fields;
  bool has_flags;
  int nullable_fields;

  std::vector<uint8_t> flags;

  std::string getFilenameX(std::string filename);

  void getFlags();

  bool skipField(const Field &field, uint8_t &ifield_for_flag_test);

  GDBTable(std::string filename);
};



class MasterTable : public GDBTable {
 public:
  std::vector< std::pair<std::string, int> > rasters;

  MasterTable(std::string filename);
};

class RasterBase : public GDBTable, public GDALGeorefPamDataset {
 private:
  std::string bandTypeToDataTypeString(std::vector< uint8_t > &band_types) const;
  std::string bandTypeToCompressionTypeString(std::vector<uint8_t> &band_types) const;
 public:
  int32_t     block_width;              //Pixel width of tiles
  int32_t     block_height;             //Pixel height of tiles
  int32_t     band_width;               //Pixel width of the band
  int32_t     band_height;              //Pixel height of the band
  double      eminx;                    //Minimum X coordinate
  double      eminy;                    //Minimum Y coordinate
  double      emaxx;                    //Maximum X coordinate
  double      emaxy;                    //Maximum Y coordinate
  double      block_origin_x;
  double      block_origin_y;
  std::string data_type;
  std::string compression_type;
  std::string name;                     //Name of the band
  int32_t     cdate;                    //Creation date
  int32_t     mdate;                    //Last modification date
  std::vector< uint8_t > band_types;
  std::array<double,6> geotransform;

  RasterBase(std::string filename);
};

class RasterProjection : public GDBTable {
 public:
  RasterProjection(std::string filename);
};

class Raster;

template<class T>
class RasterBand : public GDBTable, public GDALPamRasterBand {
 public:
  std::vector<T> geodata;

  int minpx = 0;
  int minpy = 0;
  int maxpx = 0;
  int maxpy = 0;

  RasterBand(std::string filename, const RasterBase &rb);

  void getDimensionsFromData(std::string filename, const RasterBase &rb);

  T   no_data;
  std::array<double,6> geotransform;
  std::string projection;
  int width;
  int height;
  void resize(int64_t width, int64_t height, T no_data_val);
  bool in_raster(int x, int y) const;
  T& operator()(int64_t x, int64_t y);
  T  operator()(int64_t x, int64_t y) const;
  void setAll(T val);

  GDALDataType myGDALType() const {
    if(typeid(T)==typeid(uint8_t))
      return GDT_Byte;
    else if(typeid(T)==typeid(int8_t))
      return GDT_Byte;
    else if(typeid(T)==typeid(uint16_t))
      return GDT_UInt16;
    else if(typeid(T)==typeid(int16_t))
      return GDT_Int16;
    else if(typeid(T)==typeid(uint32_t))
      return GDT_UInt32;
    else if(typeid(T)==typeid(int32_t))
      return GDT_Int32;
    else if(typeid(T)==typeid(float))
      return GDT_Float32;
    else if(typeid(T)==typeid(double))
      return GDT_Float64;
    else {
      std::cerr<<"Could not map native type '"<<typeid(T).name()<<"' to GDAL type! (Use `c++filt -t` to decode.)"<<std::endl;
      throw std::runtime_error("Could not map native data type to GDAL type!");
    }
    return GDT_Unknown;
  }

  void save(std::string filename, std::string metadata, bool compress) {
    GDALAllRegister();

    char **papszOptions = NULL;
    if(compress){
      papszOptions = CSLSetNameValue( papszOptions, "COMPRESS", "DEFLATE" );
      papszOptions = CSLSetNameValue( papszOptions, "ZLEVEL",   "6" );
    }

    if(typeid(T)==typeid(int8_t)){
      papszOptions = CSLSetNameValue( papszOptions, "PIXELTYPE", "SIGNEDBYTE" );
    }

    GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
    if(poDriver==NULL){
      std::cerr<<"Could not open GDAL driver!"<<std::endl;
      throw std::runtime_error("Could not open GDAL driver!");
    }
    GDALDataset *fout = poDriver->Create(filename.c_str(), width, height, 1, myGDALType(), papszOptions);
    if(fout==NULL){
      std::cerr<<"Could not open file '"<<filename<<"' for GDAL save!"<<std::endl;
      throw std::runtime_error("Could not open file for GDAL save!");
    }

    GDALRasterBand *oband = fout->GetRasterBand(1);
    oband->SetNoDataValue(no_data);

    //This could be used to copy metadata
    //poDstDS->SetMetadata( poSrcDS->GetMetadata() );

    //TIFFTAG_SOFTWARE
    //TIFFTAG_ARTIST
    {
      std::time_t the_time = std::time(nullptr);
      char time_str[64];
      std::strftime(time_str, sizeof(time_str), "%Y-%m-%d %H:%M:%S UTC", std::gmtime(&the_time));
      fout->SetMetadataItem("TIFFTAG_DATETIME",   time_str);
      fout->SetMetadataItem("TIFFTAG_SOFTWARE", program_identifier.c_str());

      auto out_processing_history = std::string(time_str) + " | " + program_identifier + " | ";
      if(!metadata.empty())
        out_processing_history += metadata;
      else
        out_processing_history += "Unspecified Operation";

      fout->SetMetadataItem("PROCESSING_HISTORY", out_processing_history.c_str());
    }

    //The geotransform maps each grid cell to a point in an affine-transformed
    //projection of the actual terrain. The geostransform is specified as follows:
    //    Xgeo = GT(0) + Xpixel*GT(1) + Yline*GT(2)
    //    Ygeo = GT(3) + Xpixel*GT(4) + Yline*GT(5)
    //In case of north up images, the GT(2) and GT(4) coefficients are zero, and
    //the GT(1) is pixel width, and GT(5) is pixel height. The (GT(0),GT(3))
    //position is the top left corner of the top left pixel of the raster.

    if(!geotransform.empty()){
      if(geotransform.size()!=6){
        std::cerr<<"Geotransform of output is not the right size. Found "<<geotransform.size()<<" expected 6."<<std::endl;
        throw std::runtime_error("saveGDAL(): Invalid output geotransform.");
      }

      fout->SetGeoTransform(geotransform.data());
    }

    if(!projection.empty())
      fout->SetProjection(projection.c_str());

    #ifdef DEBUG
      std::cerr<<"Filename: "<<std::setw(20)<<filename<<" Geotrans0: "<<std::setw(10)<<std::setprecision(10)<<std::fixed<<geotransform[0]<<" Geotrans3: "<<std::setw(10)<<std::setprecision(10)<<std::fixed<<geotransform[3]<< std::endl;
    #endif

    auto temp = oband->RasterIO(GF_Write, 0, 0, width, height, geodata.data(), width, height, myGDALType(), 0, 0);
    if(temp!=CE_None)
      std::cerr<<"Error writing file! Continuing in the hopes that some work can be salvaged."<<std::endl;

    GDALClose(fout);
  }
};



////////////////////////////////////////////////////////////
//UTILITY FUNCTIONS FOR READING AND MANIPULATING BINARY DATA


void bitsetToString(const std::vector< uint8_t > &bs){
  for(const auto &v: bs){
    std::bitset<8> t(v);
    std::cerr<<t<<" ";
  }
  std::cerr<<std::endl;
}

std::string hexify(int raster_num){
  std::stringstream ss;
  ss << "a" << std::setfill('0') << std::setw(8) << std::hex << raster_num << ".gdbtable";
  return ss.str();
}

template<class T>
T ReadThing(std::ifstream &fin){
  T v;
  fin.read( reinterpret_cast <char*> (&v), sizeof( T ) );
  return v;
}

uint8_t ReadByte(std::ifstream &fin){
  return ReadThing<uint8_t>(fin);
}

std::vector<uint8_t> ReadBytes(std::ifstream &fin, int count){
  std::vector<uint8_t> ret(count);
  fin.read( reinterpret_cast <char*>(ret.data()), sizeof( uint8_t )*count );
  return ret;
}

std::string ReadBytesAsString(std::ifstream &fin, int count){
  auto v = ReadBytes(fin,count);
  return std::string(v.begin(),v.end());
}

int16_t ReadInt16(std::ifstream &fin){
  #if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
    return ReadThing<int16_t>(fin);
  #else
    #pragma message "Big-endian unimplemented"
  #endif
}

int32_t ReadInt32(std::ifstream &fin){
  #if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
    return ReadThing<int32_t>(fin);
  #else
    #pragma message "Big-endian unimplemented"
  #endif
}

uint64_t ReadIndex40(std::ifstream& fin){
  #if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
  uint64_t value = ReadThing<uint32_t>(fin);
  uint64_t byte5 = ReadThing<uint8_t>(fin);
  value |= (byte5 << 32);
  return value;
  #else
  #pragma message "Big-endian unimplemented"
  #endif
}

float ReadFloat32(std::ifstream &fin){
  #if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
    return ReadThing<float>(fin);
  #else
    #pragma message "Big-endian unimplemented"
  #endif
}

double ReadFloat64(std::ifstream &fin){
  #if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
    return ReadThing<double>(fin);
  #else
    #pragma message "Big-endian unimplemented"
  #endif
}

uint64_t ReadVarUint(std::ifstream &fin){
  uint64_t shift = 0;
  uint64_t ret   = 0;
  while(true){
    uint8_t b = ReadByte(fin);
    ret      |= ((b & 0x7F) << shift);
    if( (b & 0x80)==0)
      break;
    shift += 7;
  }
  return ret;
}

void AdvanceBytes(std::ifstream &fin, int64_t count){
  fin.seekg(count, std::ios_base::cur);
}

void GotoPosition(std::ifstream &fin, int64_t pos){
  fin.seekg(pos);
}

std::string GetString(std::ifstream &fin, int nbcar=-1){
  std::string temp;

  if(nbcar==-1)
    nbcar = ReadByte(fin);

  for(int j=0;j<nbcar;j++){
    temp += ReadByte(fin);
    AdvanceBytes(fin,1);
  }
  return temp;
}

int32_t GetCount(std::ifstream &fin){
  int32_t count = ReadByte(fin);
  count += ((int)ReadByte(fin))*256;
  return count;
}


//////////////////////////////////////////////////////////
//Deal with compressed data

void Zinflate(std::vector<uint8_t> &src, std::vector<uint8_t> &dst) {
  z_stream strm  = {0};
  strm.total_in  = strm.avail_in  = src.size();
  strm.total_out = strm.avail_out = dst.size();
  strm.next_in   = (Bytef *) src.data();
  strm.next_out  = (Bytef *) dst.data();

  strm.zalloc = Z_NULL;
  strm.zfree  = Z_NULL;
  strm.opaque = Z_NULL;

  int err = -1;
  int ret = -1;

  err = inflateInit2(&strm, (15 + 32)); //15 window bits, and the +32 tells zlib to to detect if using gzip or zlib
  if(err!=Z_OK){
    inflateEnd(&strm);
    throw std::runtime_error("zlib inflateInit2 error: "+std::to_string(err));
  }

  err = inflate(&strm, Z_FINISH);
  if (err!=Z_STREAM_END) {
   inflateEnd(&strm);
   throw std::runtime_error("zlib inflate error: "+std::to_string(err));
  }

  ret = strm.total_out;
  inflateEnd(&strm);
  dst.resize(ret);
}



void Field::print() const {
  std::cout<<"Name     = "<<name      <<"\n"
           <<"Alias    = "<<alias     <<"\n"
           <<"Type     = "<<(int)type <<"\n"
           <<"Nullable = "<<nullable  <<"\n";
}


//NOTE: The following assumes that data is always stored in big endian order.
//      This should be confirmed.
template<class T>
std::vector<T> Unpack(std::vector<uint8_t> &packed, const int block_width, const int block_height){
  std::vector<T> output(block_width*block_height);
  packed.resize(sizeof(T)*block_width*block_height);

  if(std::is_same<T,float>::value){
    #if __FLOAT_WORD_ORDER__ == __ORDER_LITTLE_ENDIAN__
      for(unsigned int i=0;i<packed.size();i+=4){
        std::swap(packed[i+0],packed[i+3]);
        std::swap(packed[i+1],packed[i+2]);
      }
    #endif
  } else if(std::is_same<T,double>::value){
    #if __FLOAT_WORD_ORDER__ == __ORDER_LITTLE_ENDIAN__
      for(unsigned int i=0;i<packed.size();i+=8){
        std::swap(packed[i+0],packed[i+7]);
        std::swap(packed[i+1],packed[i+6]);
        std::swap(packed[i+2],packed[i+5]);
        std::swap(packed[i+3],packed[i+4]);
      }
    #endif
  } else if(std::is_same<T,int16_t>::value){
    #if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
      for(unsigned int i=0;i<packed.size();i+=2)
        std::swap(packed[i+0],packed[i+1]);
    #endif
  } else if(std::is_same<T,int32_t>::value){
    #if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
      for(unsigned int i=0;i<packed.size();i+=4){
        std::swap(packed[i+0],packed[i+3]);
        std::swap(packed[i+1],packed[i+2]);
      }
    #endif
  } else if(std::is_same<T,uint16_t>::value){
    #if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
      for(unsigned int i=0;i<packed.size();i+=2)
        std::swap(packed[i+0],packed[i+1]);
    #endif
  } else if(std::is_same<T,uint32_t>::value){
    #if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
      for(unsigned int i=0;i<packed.size();i+=4){
        std::swap(packed[i+0],packed[i+3]);
        std::swap(packed[i+1],packed[i+2]);
      }
    #endif
  } else if(std::is_same<T, uint8_t>::value || std::is_same<T, int8_t>::value){
    //No special unpacking needs to be done for these
  } else {
    std::cerr<<"Unimplemented type conversion for '"<<typeid(T).name()<<"'! (Use c++filt -t to decode.)"<<std::endl;
    throw std::runtime_error("Unimplemented type conversion!");
  }

  memcpy(output.data(), packed.data(), sizeof(T)*block_width*block_height);

  return output;
}


std::string GDBTable::getFilenameX(std::string filename){
  return filename.substr(0,filename.size()-1)+"x";
}

void GDBTable::getFlags(){
  if(has_flags){
    auto nremainingflags = nullable_fields;
    while(nremainingflags>0){
      auto tempflag = ReadByte(gdbtable);
      flags.push_back(tempflag);
      nremainingflags -= 8;
    }
  }
}

bool GDBTable::skipField(const Field &field, uint8_t &ifield_for_flag_test){
  if(has_flags && field.nullable){
    uint8_t test = (flags[ifield_for_flag_test >> 3] & (1 << (ifield_for_flag_test % 8)));
    ifield_for_flag_test++;
    return test!=0;
  }
  return false;
}

GDBTable::GDBTable(std::string filename){
  std::string filenamex = getFilenameX(filename);
  gdbtablx.open(filenamex, std::ios_base::in | std::ios_base::binary);

  if(!gdbtablx.good()){
    std::cerr<<"Could not find '"<<filenamex<<"'!"<<std::endl;
    throw std::runtime_error("Could not find '"+filenamex+"'!");
  }

  gdbtablx.seekg(8);
  nfeaturesx         = ReadInt32(gdbtablx);
  size_tablx_offsets = ReadInt32(gdbtablx);

  gdbtable.open(filename, std::ios_base::in | std::ios_base::binary);

  if(!gdbtable.good()){
    std::cerr<<"Could not find '"<<filename<<"'!"<<std::endl;
    throw std::runtime_error("Could not find '"+filename+"'!");
  }

  gdbtable.seekg(4);
  auto nfeatures = ReadInt32(gdbtable);

  gdbtable.seekg(32);
  auto header_offset = ReadInt32(gdbtable);

  gdbtable.seekg(header_offset);
  auto header_length = ReadInt32(gdbtable);

  AdvanceBytes(gdbtable,4);

  //1 = point
  //2 = multipoint
  //3 = polyline
  //4 = polygon
  //9 = multipatch
  auto layer_geom_type = ReadByte(gdbtable);

  AdvanceBytes(gdbtable,3);

  auto nfields = GetCount(gdbtable);

  // std::cout<<"nfeaturesx         = "<<nfeaturesx           <<"\n";
  // std::cout<<"size_tablx_offsets = "<<size_tablx_offsets   <<"\n";
  // std::cout<<"nfeatures          = "<<nfeatures            <<"\n";
  // std::cout<<"header_offset      = "<<header_offset        <<"\n";
  // std::cout<<"header_length      = "<<header_length        <<"\n";
  // std::cout<<"layer_geom_type    = "<<(int)layer_geom_type <<"\n";
  // std::cout<<"nfields            = "<<nfields              <<"\n";

  has_flags       = false;
  nullable_fields = 0;

  for(int i=0;i<nfields;i++){
    int8_t nbcar;

    auto field = Field();

    field.name     = GetString(gdbtable);
    field.alias    = GetString(gdbtable);
    field.type     = ReadByte(gdbtable);
    field.nullable = true;
    //print('type = %d (%s)' % (type, field_type_to_str(type))) //TODO

    if(field.type==6){        //ObjectID
      auto magic_byte1 = ReadByte(gdbtable);
      auto magic_byte2 = ReadByte(gdbtable);
      field.nullable   = false;
    } else if(field.type==7){ //Shape (NOTE: This is not like a shapefile, it's more like a single outline.)
      const auto magic_byte1 = ReadByte(gdbtable); //0
      const auto flag        = ReadByte(gdbtable);       //6 or 7

      if( (flag & 1)==0 )
        field.nullable = false;

      const auto wkt_len = GetCount(gdbtable);
      field.shape.wkt    = GetString(gdbtable, wkt_len/2);

      const auto magic_byte3 = ReadByte(gdbtable);

      field.shape.has_m = false;
      field.shape.has_z = false;
      if(magic_byte3==5)
        field.shape.has_z = true;
      if(magic_byte3==7){
        field.shape.has_m = true;
        field.shape.has_z = true;
      }

      field.shape.xorig   = ReadFloat64(gdbtable);
      field.shape.yorig   = ReadFloat64(gdbtable);
      field.shape.xyscale = ReadFloat64(gdbtable);
      if(field.shape.has_m){
        field.shape.morig  = ReadFloat64(gdbtable);
        field.shape.mscale = ReadFloat64(gdbtable);
      }
      if(field.shape.has_z){
        field.shape.zorig  = ReadFloat64(gdbtable);
        field.shape.zscale = ReadFloat64(gdbtable);
      }
      field.shape.xytolerance = ReadFloat64(gdbtable);
      if(field.shape.has_m)
        field.shape.mtolerance = ReadFloat64(gdbtable);
      if(field.shape.has_z)
        field.shape.ztolerance = ReadFloat64(gdbtable);

      field.shape.xmin = ReadFloat64(gdbtable);
      field.shape.ymin = ReadFloat64(gdbtable);
      field.shape.xmax = ReadFloat64(gdbtable);
      field.shape.ymax = ReadFloat64(gdbtable);

      //TODO: What is this doing?
      while(true){
        auto read5 = ReadBytes(gdbtable,5);
        if(read5[0]!=0 || (read5[1]!=1 && read5[1]!=2 && read5[1]!=3) || read5[2]!=0 || read5[3]!=0 || read5[4]!=0){
          AdvanceBytes(gdbtable,-5);
          auto datum = ReadFloat64(gdbtable);
        } else {
          for(int i=0;i<read5[1];i++)
            auto datum = ReadFloat64(gdbtable);
          break;
        }
      }

    } else if(field.type==4){ //String
      const auto width = ReadInt32(gdbtable);
      const auto flag  = ReadByte(gdbtable);
      if( (flag&1)==0 )
        field.nullable = false;
      const auto default_value_length = ReadVarUint(gdbtable);
      if( (flag&4)!=0 && default_value_length>0){
        //auto default_value = ;
        AdvanceBytes(gdbtable, default_value_length);
      }

    } else if(field.type==8){ //TODO: What is this?
      AdvanceBytes(gdbtable,1);
      const auto flag = ReadByte(gdbtable);
      if( (flag&1)==0 )
        field.nullable = false;

    } else if(field.type==9) { //Raster
      AdvanceBytes(gdbtable,1);
      const auto flag = ReadByte(gdbtable);
      if( (flag & 1)==0 )
        field.nullable = false;

      field.raster.raster_column = GetString(gdbtable);

      const auto wkt_len = GetCount(gdbtable);
      field.raster.wkt   = GetString(gdbtable,wkt_len/2);
        std::cerr<<"field.raster.wkt: "<<field.raster.wkt<<std::endl;
      #endif

      //f.read(82) //TODO: Was like this in source.

      const auto magic_byte3 = ReadByte(gdbtable);

      if(magic_byte3>0){
        field.raster.raster_has_m = false;
        field.raster.raster_has_z = false;

        if(magic_byte3==5){
          field.raster.raster_has_z = true;
        } else if(magic_byte3==7){
          field.raster.raster_has_m = true;
          field.raster.raster_has_z = true;
        }

        field.raster.raster_xorig   = ReadFloat64(gdbtable);
        field.raster.raster_yorig   = ReadFloat64(gdbtable);
        field.raster.raster_xyscale = ReadFloat64(gdbtable);

        if(field.raster.raster_has_m){
          field.raster.raster_morig  = ReadFloat64(gdbtable);
          field.raster.raster_mscale = ReadFloat64(gdbtable);
        }

        if(field.raster.raster_has_z){
          field.raster.raster_zorig  = ReadFloat64(gdbtable);
          field.raster.raster_zscale = ReadFloat64(gdbtable);
        }

        field.raster.raster_xytolerance  = ReadFloat64(gdbtable);
        if(field.raster.raster_has_m)
          field.raster.raster_mtolerance = ReadFloat64(gdbtable);
        if(field.raster.raster_has_z)
          field.raster.raster_ztolerance = ReadFloat64(gdbtable);
      }

      AdvanceBytes(gdbtable,1);

    } else if(field.type==11 || field.type==10 || field.type==12){ //UUID or XML
      const auto width = ReadByte(gdbtable);
      const auto flag  = ReadByte(gdbtable);
      if( (flag&1)==0 )
        field.nullable = false;
    } else {
      const auto width    = ReadByte(gdbtable);
      const auto flag     = ReadByte(gdbtable);
      if( (flag&1)==0 )
        field.nullable = false;

      const auto default_value_length = ReadByte(gdbtable);

      //TODO: What is this?
      if( (flag&4)!=0 ){
        if(field.type==0 && default_value_length==2)
          auto default_value = ReadInt16(gdbtable);
        else if(field.type==1 && default_value_length==4)
          auto default_value = ReadInt32(gdbtable);
        else if(field.type==2 && default_value_length==4)
          auto default_value = ReadFloat32(gdbtable);
        else if(field.type==3 && default_value_length==8)
          auto default_value = ReadFloat64(gdbtable);
        else if(field.type==5 && default_value_length==8)
          auto default_value = ReadFloat64(gdbtable);
        else
          AdvanceBytes(gdbtable, default_value_length);
      }
    }

    if(field.nullable){
      has_flags        = true;
      nullable_fields += 1;
    }

    if(field.type!=6)
      fields.push_back(field);

    //std::cout<<"\n\nField Number = "<<(fields.size()-1)<<"\n";
    //field.print();
  }
}


MasterTable::MasterTable(std::string foldername) : GDBTable(foldername + "a00000001.gdbtable") {

std::string masterTable = foldername + "a00000001.gdbtable";
  for(int f=0;f<nfeaturesx;f++){
    GotoPosition(gdbtablx, 16 + f * size_tablx_offsets);
    const auto feature_offset = ReadIndex40(gdbtablx);

    if(feature_offset==0)
      continue;

    GotoPosition(gdbtable, feature_offset);

    const auto blob_len = ReadInt32(gdbtable);

    getFlags();

    uint8_t ifield_for_flag_test = 0;
    for(unsigned int fi=0;fi<fields.size();fi++){
      if(skipField(fields[fi],ifield_for_flag_test))
        continue;

      if(fields[fi].type==1){
        const auto val = ReadInt32(gdbtable);
      } else if(fields[fi].type == 10 || fields[fi].type == 11){ //10=DatasetGUID
        const auto val = ReadBytes(gdbtable, 16);
      } else if(fields[fi].type == 4 || fields[fi].type == 12){  //String
        const auto length = ReadVarUint(gdbtable);
        const auto val    = ReadBytesAsString(gdbtable, length);
        const auto loc    = val.find("fras_ras_");
        if(loc!=std::string::npos)
          rasters.emplace_back(val.substr(9), f);
      }
    }
  }
}

RasterBase MasterTable::getRasterBase( int raster_num ){
	RasterBase rb( foldername + hexify( raster_num + 4 ) );
	return rb;
}

 vector<std::string> MasterTable::listRasters(){
  vector<std::string> rasters;
    if( rasters.size() == 0 ){
	  return bands;
    }

    for(unsigned int r  = 0; r < rasters.size(); r++ ){
	  bands.insert( rasters[r].first );
      //std::cout<<" ("<<std::hex<<std::setw(3)<<mt.rasters[r].second<<std::dec<<")"; //Latter part displays file identifier
    }
	return rasters;
  } 
  
 int MasterTable::getRasterId( std::string name ){
    if( rasters.size() == 0 ){
      return -1;
    }

    unsigned int raster_num = (unsigned int)-1;

    try {
      raster_num = std::stoi(name);
    } catch (...) {
      for( unsigned int i=0; i < rasters.size(); i++ ){
        std::locale loc;
        std::string lower;
        for( auto &s : rasters[i].first )
          s = std::tolower( s, loc );

        if( rasters[i].first == name || lower == name ){
          raster_num = i;
          break;
        }
      }
    }
    if( raster_num >= rasters.size() ){ //Note: Don't need <0 check because raster_num is unsigned
      return -1;
    }
	return raster_num;
  }
  


/* Plus 0
FID = 1
feature_offset = 498
blob_len = 96
flags = [0, 0, 252]
Type: 1
Field  sequence_nbr : 1
Type: 1
Field     raster_id : 1
Type: 4
Field          name : ""
Type: 1
Field    band_flags : 2048
Type: 1
Field    band_width : 1
Type: 1
Field   band_height : 1
Type: 1
Field    band_types : 4195328
Type: 1
Field   block_width : 128
Type: 1
Field  block_height : 128
Type: 3
Field block_origin_x : -178.500000
Type: 3
Field block_origin_y : 88.500000
Type: 3
Field         eminx : -180.000000
Type: 3
Field         eminy : 87.000000
Type: 3
Field         emaxx : -177.000000
Type: 3
Field         emaxy : 90.000000
Type: 1
Field         cdate : 1471710590
Type: 1
Field         mdate : 1471710590
Type: 1
Field          srid : 0
set([1, 3, 4, 6])
*/
RasterBase::RasterBase(std::string filename) : GDBTable(filename) {

  for(int f=0;f<nfeaturesx;f++){
    GotoPosition(gdbtablx, 16 + f * size_tablx_offsets);
    const auto feature_offset = ReadIndex40(gdbtablx);

    if(feature_offset==0)
      continue;

    GotoPosition(gdbtable, feature_offset);

    const auto blob_len = ReadInt32(gdbtable);

    getFlags();

    uint8_t ifield_for_flag_test = 0;
    for(unsigned int fi=0;fi<fields.size();fi++){
      if(skipField(fields[fi], ifield_for_flag_test))
        continue;

      if(fields[fi].type==1 && fields[fi].name=="band_types"){
        //auto val = ReadInt32(gdbtable);
        //std::cerr<<"band_types = "<<val<<std::endl;
        band_types = ReadBytes(gdbtable, 4);

        data_type        = bandTypeToDataTypeString(band_types);
        compression_type = bandTypeToCompressionTypeString(band_types);
      } else if (fields[fi].type==1) {
        const auto val = ReadInt32(gdbtable);
        if(fields[fi].name=="block_width")
          block_width = val;
        else if(fields[fi].name=="block_height")
          block_height = val;
        else if(fields[fi].name=="band_width"){
          band_width = val;
        } else if(fields[fi].name=="band_height"){
          band_height = val;
        } else if(fields[fi].name=="cdate"){
          cdate = val;
        } else if(fields[fi].name=="mdate"){
          mdate = val;
        } else {
          std::cerr<<"Skipping field '"<<fields[fi].name<<"'"<<std::endl;
        }

      } else if(fields[fi].type == 4 || fields[fi].type == 12){ //4 is a string, 12 is a UUID/XML
        const auto length = ReadVarUint(gdbtable);
        const auto val    = ReadBytesAsString(gdbtable, length);
        if(fields[fi].name=="name"){
          name = val;
        }

      } else if(fields[fi].type==3){
        const auto val = ReadFloat64(gdbtable);
        if(fields[fi].name=="block_origin_x")
          block_origin_x=val;
        else if(fields[fi].name=="block_origin_y")
          block_origin_y=val;
        else if(fields[fi].name=="eminx")
          eminx=val;
        else if(fields[fi].name=="emaxx")
          emaxx=val;
        else if(fields[fi].name=="eminy")
          eminy=val;
        else if(fields[fi].name=="emaxy")
          emaxy=val;

      }
    }
  }


  //Try to construct a geotransform. These have the form
  //  Xgeo = GT(0) + Xpixel*GT(1) + Yline*GT(2)
  //  Ygeo = GT(3) + Xpixel*GT(4) + Yline*GT(5)
  //where, in case of north up images, the GT(2) and GT(4) coefficients are
  //zero, and the GT(1) is pixel width, and GT(5) is pixel height. The
  //(GT(0),GT(3)) position is the top left corner of the top left pixel of the
  //raster.

  //TODO: This is not guaranteed to be correct
  const double pw = (emaxx-eminx)/(double)(band_width-1);
  const double ph = (emaxx-eminx)/(double)(band_width-1);
  geotransform[0] = eminx+(pw*0.5);
  geotransform[1] = pw;
  geotransform[2] = 0;
  geotransform[3] = emaxy+(ph*0.5);
  geotransform[4] = 0;
  geotransform[5] = -ph; //Arc really doesn't seem to like rasters with this value positive.

  std::cerr<<"Using geotransform (this is experimental): ";
  for(const auto &x: geotransform)
    std::cerr<<std::fixed<<std::setprecision(20)<<x<<" ";
  std::cerr<<std::endl;
}

/*
1_bit           0 4 8  0 00000000 00000100 00001000 00000000
4_bit           0 4 20 0 00000000 00000100 00100000 00000000
8_bit_signed    0 4 41 0 00000000 00000100 01000001 00000000
8_bit_unsigned  0 4 40 0 00000000 00000100 01000000 00000000
16_bit_signed   0 4 81 0 00000000 00000100 10000001 00000000
16_bit_unsigned 0 4 80 0 00000000 00000100 10000000 00000000
32_bit_signed   0 4 1  1 00000000 00000100 00000001 00000001
32_bit_float    0 4 2  1 00000000 00000100 00000010 00000001
32_bit_unsigned 0 4 0  1 00000000 00000100 00000000 00000001
64_bit          0 4 0  2 00000000 00000100 00000000 00000010
*/
std::string RasterBase::bandTypeToDataTypeString(std::vector<uint8_t> &band_types) const {
  if(band_types[2]==0x08 && band_types[3]==0x00) //00000000 00000100 00001000 00000000
    return "1bit";
  if(band_types[2]==0x20 && band_types[3]==0x00) //00000000 00000100 00100000 00000000
    return "4bit";
  if(band_types[2]==0x41 && band_types[3]==0x00) //00000000 00000100 01000001 00000000
    return "int8_t";
  if(band_types[2]==0x40 && band_types[3]==0x00) //00000000 00000100 01000000 0000000
    return "uint8_t";
  if(band_types[2]==0x81 && band_types[3]==0x00) //00000000 00000100 10000001 00000000
    return "int16_t";
  if(band_types[2]==0x80 && band_types[3]==0x00) //00000000 00000100 10000000 00000000
    return "uint16_t";
  if(band_types[2]==0x01 && band_types[3]==0x01) //00000000 00000100 00000001 00000001
    return "int32_t";
  if(band_types[2]==0x02 && band_types[3]==0x01) //00000000 00000100 00000010 00000001
    return "float32";
  if(band_types[2]==0x00 && band_types[3]==0x01) //00000000 00000100 00000000 00000001
    return "uint32_t";
  if(band_types[2]==0x00 && band_types[3]==0x02) //00000000 00000100 00000000 00000010
    return "64bit";

  std::cerr<<"Unrecognised band data type!"<<std::endl;
  throw std::runtime_error("Unrecognised band data type!");
}

std::string RasterBase::rasterType(){
  if(rb.data_type=="64bit")
    return "double";
  else if(rb.data_type=="float32")
    return "float";
  else if(rb.data_type=="uint8_t")
    return "uint8_t";
  else if(rb.data_type=="int16_t")
    return "int16_t";
  else if(rb.data_type=="int32_t")
    return "int32_t";
  else if(rb.data_type=="int8_t")
    return "int8_t";
  else if(rb.data_type=="uint16_t")
    return "uint16_t";
  else if(rb.data_type=="uint32_t")
    return "uint32_t";
  else if(rb.data_type=="4bit")
    return "uint8_t";
  else if(rb.data_type=="1bit")
    return "uint8_t";
  else
	return "unknown";
}


/*
uncompressed, float
band_types = 0 0 2  1 00000000 00000000 00000010 00000001

LZ77, float32
band_types = 0 4 2  1 00000000 00000100 00000010 00000001

JPEG, 75% quality, uint8_t
band_types = 0 8 40 0 00000000 00001000 01000000 00000000

JPEG, 23% quality, uint8_t
band_types = 0 8 40 0 00000000 00001000 01000000 00000000

JPEG2000 75% quality, int16_t
band_types = 0 c 81 0 00000000 00001100 10000001 00000000

JPEG2000 23% quality, int16_t
band_types = 0 c 81 0 00000000 00001100 10000001 00000000
*/
std::string RasterBase::bandTypeToCompressionTypeString(std::vector<uint8_t> &band_types) const {
  if(band_types[1]==0x00) //band_types = 0 0 2  1 00000000 00000000 00000010 00000001
    return "uncompressed";
  if(band_types[1]==0x04) //band_types = 0 4 2  1 00000000 00000100 00000010 00000001
    return "lz77";
  if(band_types[1]==0x08) //band_types = 0 8 40 0 00000000 00001000 01000000 00000000
    return "jpeg";
  if(band_types[1]==0x0C) //band_types = 0 c 81 0 00000000 00001100 10000001 00000000
    return "jpeg2000";

  std::cerr<<"Unrecognised band compression type!"<<std::endl;
  throw std::runtime_error("Unrecognised band compression type!");
}



RasterProjection::RasterProjection(std::string filename) : GDBTable(filename){

}


template<class T>
RasterBand<T>::RasterBand(std::string filename, const RasterBase &rb) : GDBTable(filename){
  //Determine maximum and minimum pixel coordinates from the data itself, since
  //extracting them from the metadata is not yet reliable.
  getDimensionsFromData(filename,rb);

  //NOTE: In theory, the dimensions are given by
  //    resize(rb.band_width, rb.band_height, -9999);
  //However, the pixels themselves seem to have a non-obvious coordinate scheme
  //which often takes them outside of this area. Therefore, we use
  //getDimensionsFromData() to determine the range

  ////////////////////////////////////////
  // original code here
  //resize(maxpx-minpx, maxpy-minpy, -9999);
  ////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////////////////////////
  // Original code (above) runs, but "shrink wraps" output to non-null data bounding box,
  // which is "no good" because the exported layers will not be immediately "stackable".
  // Also, there should be an option to export ALL FGDBR raster layers to a single-band GeoTIFF.
  resize(rb.band_width, rb.band_height, -9999);   
  
  //NOTE: Calculating pixel coordinates as, e.g., `col_nbr*block_width+x` may
  //yield pixel values that are not in the range [0,band_width). The following
  //offsets seemed close to working on my test data:
  //  const int xoffset = -std::abs(rb.eminx-rb.block_origin_x);
  //  const int yoffset = -std::abs(rb.emaxy-rb.block_origin_y);
  //However, there were still pixels out of range even with this offset applied
  //and the offset calculation seems sufficiently ridiculous that I don't trust
  //it. Therefore, I have set the offsets to values generated by
  //getDimensionsFromData()

    
  const int pw = (rb.emaxx-rb.eminx)/(double)rb.band_width;
  const int ph = (rb.emaxy-rb.eminy)/(double)rb.band_height;

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // This is confirmed to work for SPBC dataset (output GeoTIFF perfectly lines up with ESRI GeoTIFF).
  // Now why do we need this particular tuning factor (+13100, -10000)?
  // The values are very similar to the values of
  //   (rb.eminx-rb.block_origin_x)/pixel_width
  // and
  //   (rb.emaxy-rb.block_origin_y)/pixel_height
  // for the SPBC dataset, so there must be a way to compute the correct adjustment factor from these.
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // const int xoffset = (int)std::abs((rb.eminx-rb.block_origin_x + 13200 - 100)/pw);
  // const int yoffset = (int)std::abs((rb.emaxy-rb.block_origin_y -  9900 + 100)/ph);
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////
  // This works pefectly all 22 layers of the following large test dataset
  // https://catalogue.data.gov.bc.ca/dataset/04ad45c3-0fdc-4506-bdb4-252c45a63459
  //
  // Did extensive validation, and output lines up exactly (to pixel level) with ESRI ArcMap GeoTIFF export.
  // There might be an arithmetically simpler way to refactor this, but at least it works.
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////
  // original code
  //const int xoffset = -minpx;
  //const int yoffset = -minpy;
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////
  const int xoffset = (int)std::abs((rb.eminx-rb.block_origin_x-(int)((rb.eminx-rb.block_origin_x)/pw)-pw)/pw);
  const int yoffset = (int)std::abs((rb.emaxy-rb.block_origin_y-(int)((rb.emaxy-rb.block_origin_y)/ph)+pw)/ph);
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////  

  int skipped_points = 0;

  for(int f=0;f<nfeaturesx;f++){
    GotoPosition(gdbtablx, 16 + f * size_tablx_offsets);
    auto feature_offset = ReadIndex40(gdbtablx);

    if(feature_offset==0)
      continue;

    GotoPosition(gdbtable, feature_offset);

    auto blob_len = ReadInt32(gdbtable);

    getFlags();

    int row_nbr    = -1;
    int col_nbr    = -1;
    int rrd_factor = -1;

    //std::cerr<<"Fields.size = "<<fields.size()<<std::endl;

    uint8_t ifield_for_flag_test = 0;
    for(unsigned int fi=0;fi<fields.size();fi++){
 

      if(skipField(fields[fi], ifield_for_flag_test))
        continue;

      if(fields[fi].type==1){
        auto val = ReadInt32(gdbtable);
        //std::cout<<fields[fi].name<<" = "<<val<<std::endl;
        if(fields[fi].name=="col_nbr")
          col_nbr = val;
        else if(fields[fi].name=="row_nbr")
          row_nbr = val;
        else if(fields[fi].name=="rrd_factor")
          rrd_factor = val;
      } else if(fields[fi].type == 4 || fields[fi].type == 12){
        auto length = ReadVarUint(gdbtable);
        auto val    = ReadBytes(gdbtable, length);
      } else if(fields[fi].type==8){ //Appears to be where raster data is stored
        auto length = ReadVarUint(gdbtable);

        //Skip that which is not a base layer
        if(rrd_factor!=0){
          AdvanceBytes(gdbtable, length);
          continue;
        }

        auto val = ReadBytes(gdbtable, length);

        if(length==0)
          continue;

        std::vector<T> unpacked;

        //These two magic bytes indicate zlib compression... unless they don't,
        //since it is possible, if unlikely, for uncompressed data to begin a
        //block with these values. The `band_types` field has bits which
        //indicate compression, but appear to do so non-uniquely (lz77, lzw,
        //maybe others map to the same compression bits). Therefore, since the
        //compression indicators have not yet been entirely figured out, this
        //if- clause checks for the magic bytes and checks the length of the
        //field to determine probabilistically if compression is being used. The
        //check should be fairly robust, though, since we expect some degree of
        //compression for any non-pathological data.
        if(rb.compression_type=="lz77"){

          std::vector<uint8_t> decompressed(1000000);
          Zinflate(val, decompressed);
          decompressed.resize(sizeof(T)*rb.block_width*rb.block_height); //Drop trailer
          unpacked = Unpack<T>(decompressed, rb.block_width, rb.block_height);


        } else if(rb.compression_type=="uncompressed") {

          unpacked = Unpack<T>(val, rb.block_width, rb.block_height);
        } else {
          std::cerr<<"Unimplemented compression type!"<<std::endl;

        }



        //Save data to the numpy array
        for(int y=0;y<rb.block_height;y++)
        for(int x=0;x<rb.block_width;x++){
          int px = col_nbr*(rb.block_width )+x+xoffset;
          int py = row_nbr*(rb.block_height)+y+yoffset;
          if(in_raster(px,py)){
            operator()(px, py) = unpacked[y*(rb.block_width)+x];
          } else {
            skipped_points++;
          }
        }

      } else {
        std::cerr<<"Unrecognised field type: "<<(int)fields[fi].type<<std::endl;
      }
    }
  }


  //TODO: FGDB doesn't seem to store NoData values for the data I tested in any
  //of the gdbtables. I also found it difficult to set NoData values. At least
  //for the floating-point data I was looking at, FGDB seems to use 0xff7fffff
  //as a filler value for blocks which are not entirely filled with data. This,
  //in addition to the dataset's intrinsic NoData value need to be marked as
  //such. In the following, I assume that 0xff7fffff, 0x00, and -9999 are all
  //NoData values and translate everything to -9999 for output.

  //TODO: 0xff7fffff is how the filler value was stored in memory for the FGDB I
  //was looking at. I don't know if this is consistent across different
  //architectures or not. The following worked on my machine, but may not work
  //on yours. Some careful thinking is needed to improve this section.
  const uint32_t arc_no_data      = 0xff7fffff;
  const T *translated_arc_no_data = reinterpret_cast<const T*>(&arc_no_data);

  no_data = -9999; //TODO: This cannot always be NoData.

  for(uint64_t i=0;i<geodata.size();i++)
    if(geodata[i]==*translated_arc_no_data)
      geodata[i] = no_data;
    else if(geodata[i]==0) //TODO: Surely this is not always NoData?
      geodata[i] = no_data;
}


template<class T>
void RasterBand<T>::getDimensionsFromData(std::string filename, const RasterBase &rb){
  minpx = std::numeric_limits<int>::max();
  minpy = std::numeric_limits<int>::max();
  maxpx = std::numeric_limits<int>::min();
  maxpy = std::numeric_limits<int>::min();

  for(int f=0;f<nfeaturesx;f++){
    GotoPosition(gdbtablx, 16 + f * size_tablx_offsets);
    auto feature_offset = ReadIndex40(gdbtablx);

    if(feature_offset==0)
      continue;

    GotoPosition(gdbtable, feature_offset);

    auto blob_len = ReadInt32(gdbtable);

    getFlags();

    int row_nbr    = -1;
    int col_nbr    = -1;
    int rrd_factor = -1;

    uint8_t ifield_for_flag_test = 0;
    for(unsigned int fi=0;fi<fields.size();fi++){
      if(skipField(fields[fi], ifield_for_flag_test))
        continue;

      if(fields[fi].type==1){
        auto val = ReadInt32(gdbtable);
        if(fields[fi].name=="col_nbr")
          col_nbr = val;
        else if(fields[fi].name=="row_nbr")
          row_nbr = val;
        else if(fields[fi].name=="rrd_factor")
          rrd_factor = val;
      } else if(fields[fi].type == 4 || fields[fi].type == 12){
        const auto length = ReadVarUint(gdbtable);
        const auto val    = ReadBytes(gdbtable, length);
      } else if(fields[fi].type==8){ //Appears to be where raster data is stored
        if(rrd_factor!=0)
          continue;

        const int px = col_nbr*rb.block_width;
        const int py = row_nbr*rb.block_height;

        minpx = std::min(minpx,px);
        minpy = std::min(minpy,py);
        maxpx = std::max(maxpx,px+rb.block_width);
        maxpy = std::max(maxpy,py+rb.block_height);

        break;

      } else {
        std::cerr<<"Unrecognised field type: "<<(int)fields[fi].type<<std::endl;
      }
    }
  }
}


template<class T>
void RasterBand<T>::resize(int64_t width, int64_t height, T no_data_val){
  this->width  = width;
  this->height = height;
  std::cerr<<"Allocating "<<sizeof(T)<<"x"<<width<<"x"<<height<<" = "<<(sizeof(T)*width*height)<<std::endl;
  geodata.resize(width*height);
  std::fill(geodata.begin(), geodata.end(), no_data_val);
}

template<class T>
bool RasterBand<T>::in_raster(int x, int y) const {
  //std::cerr<<std::dec<<x<<" "<<width<<" "<<y<<" "<<height<<std::endl;
  return 0<=x && x<width && 0<=y && y<height;
}

template<class T>
T& RasterBand<T>::operator()(int64_t x, int64_t y){
  return geodata[y*width+x];
}

template<class T>
T RasterBand<T>::operator()(int64_t x, int64_t y) const {
  return geodata[y*width+x];
}

template<class T>
void RasterBand<T>::setAll(T val){
  std::fill(geodata.begin(),geodata.end(),val);
}



template<class T>
void ExportTypedRasterToGeoTIFF(std::string operation, std::string basename, int raster_num, std::string outputname){


class Raster(string sourceDb){

  GDBTable        bt;
  RasterProjection rp;
  RasterBand<T>  rd;
  RasterBase     rb;

}


void Raster::readRaster( std::string basename, int raster_num ){

  bt(basename+hexify(raster_num));
  rp(basename+hexify(raster_num+1));
  rb(basename+hexify(raster_num+4)); //Get the fras_bnd file 
  rd(basename+hexify(raster_num+3), rb);

  for(const auto &f: bt.fields)
    if(f.type==9){
      if(rd.projection!="")
        std::cerr<<"Ambiguity in which WKT projection to use. Using the last one."<<std::endl;
      rd.projection = f.raster.wkt;
    }

  rd.geotransform = rb.geotransform;
}



  rd.save(outputname, operation, false);
}

void ExportRasterToGeoTIFF(std::string operation, std::string basename, int raster_num, std::string outputname){
  RasterBase rb(basename+hexify(raster_num+4)); //Get the fras_bnd file


    
}

#include <locale>

int main(int argc, char **argv){
  std::string operation;
  if(argc!=2 && argc!=4){
    std::cerr<<"Syntax A: "<<argv[0]<<" <File Geodatabase>"<<std::endl;
    std::cerr<<"Syntax B: "<<argv[0]<<" <File Geodatabase> <Raster> <Output Name>"<<std::endl;
    std::cerr<<"\n";
    std::cerr<<"Syntax A will list all of the rasters in the data set along with selection numbers for use with Syntax B.\n";
    std::cerr<<"Syntax B, given an FGDB and raster selection number, will output the raster to the indicated output file.\n";
    std::cerr<<"Syntax B also accepts raster names, as listed by Syntax A, as inputs for <Raster>\n";
    std::cerr<<"\nNOTE: The geodatabase path must end with a slash!\n";
    std::cerr<<"EXAMPLE: ./arc_raster.exe path/to/geodatabase.gdb/ dem03 /z/out.tif\n";
    return -1;
  }

  std::string basename = argv[1];

  MasterTable mt(basename+"a00000001.gdbtable");

    ExportRasterToGeoTIFF(operation, basename, mt.rasters.at(raster_num).second, std::string(argv[3]));
  }

  return 0;
}
