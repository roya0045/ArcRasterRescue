

#include "gdal_pam.h"
#include "ogr_spatialref.h"

#include <algorithm>
#include <mutex>

#include "openGDB_headers.h"

using namespace OPENGDB_IMF_NAMESPACE;
using namespace IMATH_NAMESPACE;

extern "C" CPL_DLL void GDALRegister_GDB();

static  const char* const apszCompressions[] = {
    "NONE",
    "RLE",
    "ZIPS",
    "ZIP",
    "PIZ",
    "PXR24",
    "B44",
    "B44A",
    "DWAA",
    "DWAB",
};

#pragma once

#include <gdal_priv.h>

#include <algorithm>
#include <cassert>
#include <chrono>
#include <fstream>
#include <iostream>
#include <string>
#include <typeinfo>
#include <vector>

#ifndef GIT_HASH
  #pragma message "Compiling without a git hash!"
  const std::string git_hash = "NO HASH SPECIFIED!";
#else
  const std::string git_hash = std::string(GIT_HASH).substr(0,16);
#endif

#ifndef COMPILE_TIME
  #pragma message "Compiling without UTC compile time falling back to local!"
  const std::string compilation_datetime = __DATE__ " " __TIME__;
#else
  const std::string compilation_datetime = COMPILE_TIME;
#endif

const std::string program_url = "github.com/r-barnes/ArcRasterRescue";

///Richdem vX.X.X
const std::string program_name = "Arc Raster Rescue";

///Richdem vX.X.X (hash=GIT HASH, compiled=COMPILATION DATE TIME)
const std::string program_identifier = program_name + " (url="+program_url+", hash=" + git_hash + ", compiled="+compilation_datetime + ")";

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


class BaseTable {
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

  BaseTable(std::string filename);
};



class MasterTable : public GDALPamDataset {
 public:
  std::vector< std::pair<std::string, int> > rasters;

  MasterTable(std::string filename);
};

class RasterBase : public GDALPamRasterBand {
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

  BaseTable(std::string filename);
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

class RasterProjection : public BaseTable {
 public:
  RasterProjection(std::string filename);
};

class Raster;

template<class T>
class RasterData : public BaseTable {
 public:
  std::vector<T> geodata;

  int minpx = 0;
  int minpy = 0;
  int maxpx = 0;
  int maxpy = 0;

  RasterData(std::string filename, const RasterBase &rb);

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


void ExportRasterToGeoTIFF(std::string operation, std::string basename, int raster_num, std::string outputname);


/************************************************************************/
/*                       GDALGDBDataset()                               */
/************************************************************************/

class GDALGDBDataset final: public GDALPamDataset
{
        friend class GDALGDBRasterBand;
        friend class GDALGDBPreviewRasterBand;
        friend class GDALGDBRGBARasterBand;

        // Keep stream before others, so that it is destroyed last
        std::unique_ptr<IStream> m_pIStream{};

        std::unique_ptr<TiledInputPart> m_pTiledIP{};
        std::unique_ptr<InputPart> m_pIP{};

        std::unique_ptr<MultiPartInputFile> m_pMPIF{};
        std::unique_ptr<RgbaInputFile> m_pRGBAIF{};

        std::vector<Rgba> m_rgbaBuffer{};
        int m_nRGBABufferLine = -1;
        int m_iPart = 0;
        int m_nDWMinX = 0;
        int m_nDWMinY = 0;
        GDALGDBDataset* m_poParent = nullptr;
        int m_iLevel = 0;
        std::vector<std::unique_ptr<GDALGDBDataset>> m_apoOvrDS{};
        OGRSpatialReference m_oSRS{};
        double m_adfGT[6] = {0,1,0,0,0,1};
        bool m_bHasGT = false;

    public:
        GDALGDBDataset() = default;
        ~GDALGDBDataset();

        const OGRSpatialReference* GetSpatialRef() const override;
        CPLErr GetGeoTransform(double* adfGT) override;

        static int Identify(GDALOpenInfo* poOpenInfo);
        static GDALDataset* Open(GDALOpenInfo* poOpenInfo);
        static GDALDataset *Create( const char * pszFilename,
                                    int nXSize, int nYSize, int nBands,
                                    GDALDataType eType, char ** papszOptions );
        static GDALDataset *CreateCopy( const char * pszFilename,
                                    GDALDataset *poSrcDS,
                                    int bStrict, char ** papszOptions,
                                    GDALProgressFunc pfnProgress,
                                    void * pProgressData );
};

/************************************************************************/
/*                       GDALGDBRasterBand()                            */
/************************************************************************/

class GDALGDBRasterBand final: public GDALPamRasterBand
{
        friend class GDALGDBDataset;

        GDALColorInterp m_eInterp = GCI_Undefined;
        std::string m_osChannelName;

    protected:
        CPLErr IReadBlock(int, int, void*) override;

    public:
        GDALGDBRasterBand(GDALGDBDataset* poDSIn, int nBandIn,
                          const std::string& channelName,
                          PixelType pixelType,
                          int nBlockXSizeIn, int nBlockYSizeIn);

        GDALColorInterp GetColorInterpretation() override { return m_eInterp; }
        int GetOverviewCount() override;
        GDALRasterBand* GetOverview(int) override;
};

/************************************************************************/
/*                         GDALGDBRasterBand()                          */
/************************************************************************/

GDALGDBRasterBand::GDALGDBRasterBand(GDALGDBDataset* poDSIn, int nBandIn,
                                     const std::string& channelName,
                                     PixelType pixelType,
                                     int nBlockXSizeIn, int nBlockYSizeIn):
    m_osChannelName(channelName)
{
    poDS = poDSIn;
    nBand = nBandIn;
    nRasterXSize = poDSIn->GetRasterXSize();
    nRasterYSize = poDSIn->GetRasterYSize();
    nBlockXSize = nBlockXSizeIn;
    nBlockYSize = nBlockYSizeIn;
    eDataType = (pixelType == UINT) ? GDT_UInt32 : GDT_Float32;
}

/************************************************************************/
/*                          IReadBlock()                                */
/************************************************************************/

CPLErr GDALGDBRasterBand::IReadBlock(int nBlockXOff, int nBlockYOff,
                                     void* pImage)
{
    auto poGDS = cpl::down_cast<GDALGDBDataset*>(poDS);
    try
    {
        FrameBuffer fb;
        const size_t sizeOfElt = sizeof(float); // sizeof(uint32) as well
        const auto slice =
            Slice(eDataType == GDT_Float32 ? FLOAT : UINT,
                    static_cast<char*>(pImage) -
                    (poGDS->m_nDWMinX + nBlockXOff * nBlockXSize +
                        static_cast<size_t>(
                        poGDS->m_nDWMinY + nBlockYOff * nBlockYSize) * nBlockXSize) * sizeOfElt,
                    sizeOfElt,
                    sizeOfElt * nBlockXSize);
        fb.insert(m_osChannelName, slice);

        if( poGDS->m_pIP )
        {
            poGDS->m_pIP->setFrameBuffer(fb);
            poGDS->m_pIP->readPixels(poGDS->m_nDWMinY + nBlockYOff);
        }
        else
        {
            auto tiledIP = poGDS->m_poParent ?
                poGDS->m_poParent->m_pTiledIP.get() : poGDS->m_pTiledIP.get();
            CPLAssert( tiledIP );
            tiledIP->setFrameBuffer(fb);
            tiledIP->readTile(nBlockXOff, nBlockYOff, poGDS->m_iLevel);
        }
        return CE_None;
    }
    catch( const std::exception& e )
    {
        if( strstr(e.what(), "is missing") )
        {
            CPLDebug("GDB", "%s", e.what());
            memset(pImage, 0,
                   static_cast<size_t>(nBlockXSize) * nBlockYSize *
                   GDALGetDataTypeSizeBytes(eDataType));
            return CE_None;
        }
        CPLError(CE_Failure, CPLE_AppDefined, "OpenGDB: %s", e.what());
    }

    return CE_Failure;
}

/************************************************************************/
/*                          GetOverviewCount()                          */
/************************************************************************/

int GDALGDBRasterBand::GetOverviewCount()
{
    auto poGDS = cpl::down_cast<GDALGDBDataset*>(poDS);
    return static_cast<int>(poGDS->m_apoOvrDS.size());
}

/************************************************************************/
/*                            GetOverview()                             */
/************************************************************************/

GDALRasterBand* GDALGDBRasterBand::GetOverview(int iOvr)
{
    if( iOvr < 0 || iOvr >= GetOverviewCount() )
        return nullptr;
    auto poGDS = cpl::down_cast<GDALGDBDataset*>(poDS);
    return poGDS->m_apoOvrDS[iOvr]->GetRasterBand(nBand);
}

/************************************************************************/
/*                      GDALGDBPreviewRasterBand                        */
/************************************************************************/

class GDALGDBPreviewRasterBand final: public GDALPamRasterBand
{
        friend class GDALGDBDataset;

        std::string m_osChannelName;

    protected:
        CPLErr IReadBlock(int, int, void*) override;
        GDALColorInterp GetColorInterpretation() override {
            return static_cast<GDALColorInterp>(GCI_RedBand + nBand - 1); }

    public:
        GDALGDBPreviewRasterBand(GDALGDBDataset* poDSIn, int nBandIn);
};

/************************************************************************/
/*                      GDALGDBPreviewRasterBand()                      */
/************************************************************************/

GDALGDBPreviewRasterBand::GDALGDBPreviewRasterBand(GDALGDBDataset* poDSIn,
                                                   int nBandIn)
{
    poDS = poDSIn;
    nBand = nBandIn;
    nRasterXSize = poDSIn->GetRasterXSize();
    nRasterYSize = poDSIn->GetRasterYSize();
    nBlockXSize = nRasterXSize;
    nBlockYSize = 1;
    eDataType = GDT_Byte;
}

/************************************************************************/
/*                          IReadBlock()                                */
/************************************************************************/

CPLErr GDALGDBPreviewRasterBand::IReadBlock(int, int nBlockYOff,
                                            void* pImage)
{
    auto poGDS = cpl::down_cast<GDALGDBDataset*>(poDS);
    try
    {
        const auto& header = poGDS->m_pMPIF->header(poGDS->m_iPart);
        const auto& preview = header.previewImage();
        GDALCopyWords(reinterpret_cast<const GByte*>(preview.pixels() +
                        nBlockYOff * nRasterXSize) + nBand - 1,
                      GDT_Byte, 4,
                      pImage, GDT_Byte, 1,
                      nRasterXSize);
        return CE_None;
    }
    catch( const std::exception& e )
    {
        CPLError(CE_Failure, CPLE_AppDefined, "OpenGDB: %s", e.what());
    }

    return CE_Failure;
}

/************************************************************************/
/*                      GDALGDBRGBARasterBand                        */
/************************************************************************/

class GDALGDBRGBARasterBand final: public GDALPamRasterBand
{
        friend class GDALGDBDataset;

        std::string m_osChannelName;

    protected:
        CPLErr IReadBlock(int, int, void*) override;
        GDALColorInterp GetColorInterpretation() override {
            return static_cast<GDALColorInterp>(GCI_RedBand + nBand - 1); }

    public:
        GDALGDBRGBARasterBand(GDALGDBDataset* poDSIn, int nBandIn);
};

/************************************************************************/
/*                      GDALGDBRGBARasterBand()                      */
/************************************************************************/

GDALGDBRGBARasterBand::GDALGDBRGBARasterBand(GDALGDBDataset* poDSIn,
                                             int nBandIn)
{
    poDS = poDSIn;
    nBand = nBandIn;
    nRasterXSize = poDSIn->GetRasterXSize();
    nRasterYSize = poDSIn->GetRasterYSize();
    nBlockXSize = nRasterXSize;
    nBlockYSize = 1;
    eDataType = GDT_Float32;

}

/************************************************************************/
/*                          IReadBlock()                                */
/************************************************************************/

CPLErr GDALGDBRGBARasterBand::IReadBlock(int, int nBlockYOff, void* pImage)
{
    auto poGDS = cpl::down_cast<GDALGDBDataset*>(poDS);
    try
    {
        if( nBlockYOff != poGDS->m_nRGBABufferLine )
        {
            poGDS->m_rgbaBuffer.resize(nRasterXSize);
            poGDS->m_pRGBAIF->setFrameBuffer(
                poGDS->m_rgbaBuffer.data() - ((poGDS->m_nDWMinY + nBlockYOff) *
                        static_cast<size_t>(nRasterXSize) + poGDS->m_nDWMinX),
                1,
                nRasterXSize);
            poGDS->m_pRGBAIF->readPixels(poGDS->m_nDWMinY + nBlockYOff);
        }
        if( nBand == 1 )
        {
            for( int i = 0; i < nRasterXSize; i++ )
            {
                static_cast<float*>(pImage)[i] = poGDS->m_rgbaBuffer[i].r;
            }
        }
        else if( nBand == 2 )
        {
            for( int i = 0; i < nRasterXSize; i++ )
            {
                static_cast<float*>(pImage)[i] = poGDS->m_rgbaBuffer[i].g;
            }
        }
        else if( nBand == 3 )
        {
            for( int i = 0; i < nRasterXSize; i++ )
            {
                static_cast<float*>(pImage)[i] = poGDS->m_rgbaBuffer[i].b;
            }
        }
#ifdef unused
        else
        {
            for( int i = 0; i < nRasterXSize; i++ )
            {
                static_cast<float*>(pImage)[i] = poGDS->m_rgbaBuffer[i].a;
            }
        }
#endif
        poGDS->m_nRGBABufferLine = nBlockYOff;
        return CE_None;
    }
    catch( const std::exception& e )
    {
        if( strstr(e.what(), "is missing") )
        {
            CPLDebug("GDB", "%s", e.what());
            memset(pImage, 0,
                   static_cast<size_t>(nBlockXSize) * nBlockYSize *
                   GDALGetDataTypeSizeBytes(eDataType));
            return CE_None;
        }
        CPLError(CE_Failure, CPLE_AppDefined, "OpenGDB: %s", e.what());
    }

    return CE_Failure;
}

/************************************************************************/
/*                         ~GDALGDBDataset()                            */
/************************************************************************/

GDALGDBDataset::~GDALGDBDataset() = default;

/************************************************************************/
/*                          GetSpatialRef()                             */
/************************************************************************/

const OGRSpatialReference* GDALGDBDataset::GetSpatialRef() const
{
    const auto* poPamSRS = GDALPamDataset::GetSpatialRef();
    if( poPamSRS )
        return poPamSRS;
    return m_oSRS.IsEmpty() ? nullptr : &m_oSRS;
}

/************************************************************************/
/*                         GetGeoTransform()                            */
/************************************************************************/

CPLErr GDALGDBDataset::GetGeoTransform(double* adfGT)
{
    if( GDALPamDataset::GetGeoTransform(adfGT) == CE_None )
    {
        return CE_None;
    }
    memcpy(adfGT, m_adfGT, 6 * sizeof(double));
    return m_bHasGT ? CE_None : CE_Failure;
}

/************************************************************************/
/*                            Identify()                                */
/************************************************************************/

int GDALGDBDataset::Identify(GDALOpenInfo* poOpenInfo)
{
    if( STARTS_WITH_CI(poOpenInfo->pszFilename, "GDB:") )
        return true;

    // Check magic number
    return poOpenInfo->fpL != nullptr &&
           poOpenInfo->nHeaderBytes >= 4 &&
           poOpenInfo->pabyHeader[0] == 0x76 &&
           poOpenInfo->pabyHeader[1] == 0x2f &&
           poOpenInfo->pabyHeader[2] == 0x31 &&
           poOpenInfo->pabyHeader[3] == 0x01;
}

/************************************************************************/
/*                           GDALGDBIOStream                            */
/************************************************************************/

class GDALGDBIOStreamException final: public std::exception
{
    std::string m_msg;

    public:
        explicit GDALGDBIOStreamException(const std::string& msg): m_msg(msg) {}
        const char* what() const noexcept override { return m_msg.c_str(); }
};

#if OPENGDB_VERSION_MAJOR < 3
typedef Int64 IoInt64Type;
#else
typedef uint64_t IoInt64Type;
#endif

class GDALGDBIOStream final: public IStream, public OStream
{
  public:

    GDALGDBIOStream (VSILFILE* fp, const char* filename):
            IStream (filename), OStream (filename), m_fp (fp) {}
    ~GDALGDBIOStream() { VSIFCloseL(m_fp); }

    virtual bool        read (char c[/*n*/], int n) override;
    virtual void        write (const char c[/*n*/], int n) override;
    virtual IoInt64Type tellg () override;
    virtual IoInt64Type tellp () override { return tellg(); }
    virtual void        seekg (IoInt64Type pos) override;
    virtual void        seekp (IoInt64Type pos) override { return seekg(pos); }

  private:
    VSILFILE* m_fp;
};

bool GDALGDBIOStream::read (char c[/*n*/], int n)
{
    if( static_cast<int>(VSIFReadL(c, 1, n, m_fp)) != n ) {
        if( VSIFEofL(m_fp) )
        {
            throw GDALGDBIOStreamException(
                CPLSPrintf("Unexpected end of file. Cannot read %d bytes", n));
        }
        else
        {
            throw GDALGDBIOStreamException(
                CPLSPrintf("cannot read %d bytes", n));
        }
    }
    return VSIFEofL(m_fp) != 0;
}

void GDALGDBIOStream::write (const char c[/*n*/], int n)
{
    if( static_cast<int>(VSIFWriteL(c, 1, n, m_fp)) != n ) {
        throw GDALGDBIOStreamException(CPLSPrintf("cannot write %d bytes", n));
    }
}

IoInt64Type GDALGDBIOStream::tellg ()
{
    return static_cast<IoInt64Type>(VSIFTellL(m_fp));
}

void GDALGDBIOStream::seekg (IoInt64Type pos)
{
    VSIFSeekL(m_fp, static_cast<vsi_l_offset>(pos), SEEK_SET);
}

/************************************************************************/
/*                           setNumThreads()                            */
/************************************************************************/

static void setNumThreads()
{
    static std::mutex mutex;
    std::lock_guard<std::mutex> oLock(mutex);
    static bool bSet = false;
    if( !bSet )
    {
        bSet = true;
        setGlobalThreadCount(CPLGetNumCPUs());
    }
}

/************************************************************************/
/*                               Open()                                 */
/************************************************************************/

GDALDataset* GDALGDBDataset::Open(GDALOpenInfo* poOpenInfo)
{
    if( !Identify(poOpenInfo) )
        return nullptr;
    if( poOpenInfo->eAccess == GA_Update )
    {
        CPLError(CE_Failure, CPLE_NotSupported,
                 "Update of existing GDB file not supported");
        return nullptr;
    }

    CPLString osFilename(poOpenInfo->pszFilename);
    int iPart = 0;
    bool bIsPreview = false;
    VSILFILE* fp;
    if( STARTS_WITH_CI(poOpenInfo->pszFilename, "GDB:") )
    {
        bIsPreview = STARTS_WITH_CI(poOpenInfo->pszFilename, "GDB:PREVIEW:");
        const char* pszPartPos = bIsPreview ?
                poOpenInfo->pszFilename + strlen("GDB:PREVIEW:") :
                poOpenInfo->pszFilename + strlen("GDB:");
        const char* pszNextColumn = strchr(pszPartPos, ':');
        if( pszNextColumn == nullptr )
            return nullptr;
        iPart = atoi(pszPartPos);
        if( iPart <= 0 )
            return nullptr;
        osFilename = pszNextColumn + 1;
        fp = VSIFOpenL(osFilename, "rb");
        if( fp == nullptr )
            return nullptr;
    }
    else
    {
        fp = poOpenInfo->fpL;
        poOpenInfo->fpL = nullptr;
    }

    try
    {
        auto poDS = std::unique_ptr<GDALGDBDataset>(new GDALGDBDataset());
        poDS->m_pIStream.reset(new GDALGDBIOStream(fp, osFilename));
        poDS->m_pMPIF.reset(new MultiPartInputFile(*poDS->m_pIStream));
        if( iPart > 0 && iPart > poDS->m_pMPIF->parts() )
            return nullptr;

        if( iPart > 0 || poDS->m_pMPIF->parts() == 1 )
        {
            iPart = iPart > 0 ? iPart-1 : 0;
            poDS->m_iPart = iPart;

            const auto& header = poDS->m_pMPIF->header(iPart);
            if( bIsPreview )
            {
                if( !header.hasPreviewImage() )
                    return nullptr;
                for( int i = 1; i <= 4; i++ )
                {
                    const auto& preview = header.previewImage();
                    poDS->nRasterXSize = preview.width();
                    poDS->nRasterYSize = preview.height();
                    poDS->SetBand(i,
                                  new GDALGDBPreviewRasterBand(poDS.get(), i));
                }
                return poDS.release();
            }

            const auto& dataWindow = header.dataWindow();
            poDS->m_nDWMinX = dataWindow.min.x;
            poDS->m_nDWMinY = dataWindow.min.y;
            poDS->nRasterXSize = 1 + dataWindow.max.x - dataWindow.min.x;
            poDS->nRasterYSize = 1 + dataWindow.max.y - dataWindow.min.y;
            const auto &channels = header.channels();
            int i = 0;
            bool BGR = true;
            bool ABGR = true;
            bool BYRYY = true;
            PixelType samePixelType = NUM_PIXELTYPES;
            for (auto iter = channels.begin();
                      iter != channels.end(); ++iter, ++i)
            {
                const Channel &channel = iter.channel();
                const std::string name(iter.name());
                if( i == 0 )
                    samePixelType = channel.type;
                else if( samePixelType != channel.type )
                {
                    ABGR = false;
                    BGR = false;
                }

                if( i == 0 && name != "B" )
                    BGR = false;
                else if( i == 1 && name != "G" )
                    BGR = false;
                else if( i == 2 && name != "R" )
                    BGR = false;

                if( i == 0 && name != "A" )
                    ABGR = false;
                else if( i == 1 && name != "B" )
                    ABGR = false;
                else if( i == 2 && name != "G" )
                    ABGR = false;
                else if( i == 3 && name != "R" )
                    ABGR = false;

                if( i == 0 && name != "BY" )
                    BYRYY = false;
                else if( i == 1 && name != "RY" )
                    BYRYY = false;
                else if( i == 2 && name != "Y" )
                    BYRYY = false;
            }
            BGR &= (i == 3);
            ABGR &= (i == 4);
            BYRYY &= iPart == 0 && (i == 3);
            int nBlockXSize = poDS->nRasterXSize;
            int nBlockYSize = 1;
            if( header.hasTileDescription() )
            {
                const auto& tileDesc = header.tileDescription();
                nBlockXSize = tileDesc.xSize;
                nBlockYSize = tileDesc.ySize;
                poDS->m_pTiledIP.reset(new TiledInputPart(*poDS->m_pMPIF, iPart));
            }
            else if( BYRYY )
            {
                poDS->m_pIStream->seekg(0);
                poDS->m_pRGBAIF.reset(new RgbaInputFile(*poDS->m_pIStream));
            }
            else
            {
                poDS->m_pIP.reset(new InputPart(*poDS->m_pMPIF, iPart));
            }
            if( BYRYY )
            {
                for( i = 1; i <= 3; i++ )
                {
                    poDS->SetBand(i,
                                  new GDALGDBRGBARasterBand(poDS.get(), i));
                }
                poDS->SetMetadataItem("INTERLEAVE", "PIXEL",
                                      "IMAGE_STRUCTURE");
                poDS->SetMetadataItem("SOURCE_COLOR_SPACE", "YCbCr",
                                      "IMAGE_STRUCTURE");
            }
            else if( BGR || ABGR )
            {
                const int nBands = i;
                i = 0;
                for (auto iter = channels.begin();
                          iter != channels.end(); ++iter, ++i)
                {
                    auto poBand = new GDALGDBRasterBand(poDS.get(), nBands - i,
                                                        iter.name(),
                                                        samePixelType,
                                                        nBlockXSize,
                                                        nBlockYSize);
                    poBand->m_eInterp =
                        static_cast<GDALColorInterp>(GCI_RedBand + nBands - 1 - i);
                    poDS->SetBand(nBands - i, poBand);
                }
            }
            else
            {
                i = 0;
                for (auto iter = channels.begin();
                          iter != channels.end(); ++iter, ++i)
                {
                    const Channel &channel = iter.channel();
                    auto poBand = new GDALGDBRasterBand(poDS.get(), i+1,
                                                        iter.name(),
                                                        channel.type,
                                                        nBlockXSize,
                                                        nBlockYSize);
                    const std::string name(iter.name());
                    if( name != CPLSPrintf("Band%d", i+1) )
                        poBand->SetDescription(name.c_str());
                    if( name == "B" )
                        poBand->m_eInterp = GCI_BlueBand;
                    else if( name == "G" )
                        poBand->m_eInterp = GCI_GreenBand;
                    else if( name == "R" )
                        poBand->m_eInterp = GCI_RedBand;
                    else if( name == "A" )
                        poBand->m_eInterp = GCI_AlphaBand;
                    else if( name == "Y" )
                        poBand->m_eInterp = GCI_GrayIndex;
                    poDS->SetBand(i+1, poBand);
                }
            }

            if( poDS->m_pTiledIP && !BYRYY &&
                // Not completely clear on tiling & overviews would work
                // on dataWindow.min != 0, so exclude that for now
                dataWindow.min.x == 0 && dataWindow.min.y == 0 )
            {
                int nLevels = std::min(poDS->m_pTiledIP->numXLevels(),
                                       poDS->m_pTiledIP->numYLevels());
                for( int iLevel = 1; iLevel < nLevels; iLevel++ )
                {
                    const int nOvrWidth = poDS->m_pTiledIP->levelWidth(iLevel);
                    const int nOvrHeight = poDS->m_pTiledIP->levelHeight(iLevel);
                    if( nOvrWidth < 128 && nOvrHeight < 128 )
                    {
                        break;
                    }
                    auto poOvrDS =
                        std::unique_ptr<GDALGDBDataset>(new GDALGDBDataset());
                    // coverity[escape]
                    poOvrDS->m_poParent = poDS.get();
                    poOvrDS->m_iLevel = iLevel;
                    poOvrDS->nRasterXSize = nOvrWidth;
                    poOvrDS->nRasterYSize = nOvrHeight;
                    poDS->m_apoOvrDS.push_back(std::move(poOvrDS));
                    i = 0;
                    for (auto iter = channels.begin();
                              iter != channels.end(); ++iter, ++i)
                    {
                        const Channel &channel = iter.channel();
                        auto poBand = new GDALGDBRasterBand(
                            poDS->m_apoOvrDS.back().get(), i+1,
                            iter.name(),
                            channel.type,
                            nBlockXSize,
                            nBlockYSize);
                        poDS->m_apoOvrDS.back()->SetBand(i+1, poBand);
                    }
                }
            }

            for (auto iter = header.begin(); iter != header.end(); ++iter)
            {
                const Attribute *attr = &iter.attribute();
                const StringAttribute *stringAttr =
                                dynamic_cast <const StringAttribute *>(attr);
                const M33dAttribute* m33DAttr =
                                dynamic_cast <const M33dAttribute *>(attr);
                if ( stringAttr && strcmp(iter.name(), "gdal:crsWkt") == 0)
                {
                    poDS->m_oSRS.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
                    poDS->m_oSRS.importFromWkt(stringAttr->value().c_str());
                }
                else if ( m33DAttr && strcmp(iter.name(), "gdal:geoTransform") == 0)
                {
                    poDS->m_bHasGT = true;
                    poDS->m_adfGT[0] = m33DAttr->value()[0][2];
                    poDS->m_adfGT[1] = m33DAttr->value()[0][0];
                    poDS->m_adfGT[2] = m33DAttr->value()[0][1];
                    poDS->m_adfGT[3] = m33DAttr->value()[1][2];
                    poDS->m_adfGT[4] = m33DAttr->value()[1][0];
                    poDS->m_adfGT[5] = m33DAttr->value()[1][1];
                }
                else if ( stringAttr && STARTS_WITH(iter.name(), "gdal:") )
                {
                    poDS->SetMetadataItem(iter.name() + strlen("gdal:"),
                                          stringAttr->value().c_str());
                }
                else if ( stringAttr && strcmp(iter.name(), "type") != 0)
                {
                    poDS->SetMetadataItem(iter.name(),
                                          stringAttr->value().c_str());
                }
            }

            const auto& compression = header.compression();
            if( compression == NO_COMPRESSION )
            {
                // nothing
            }
            else if( compression < CPL_ARRAYSIZE(apszCompressions) )
            {
                poDS->SetMetadataItem("COMPRESSION",
                                      apszCompressions[compression],
                                      "IMAGE_STRUCTURE");
            }
            else
            {
                CPLDebug("GDB", "Unknown compression method: %d", compression);
            }

            if( header.hasPreviewImage() )
            {
                CPLStringList aosSubDS;
                aosSubDS.SetNameValue("SUBDATASET_1_NAME",
                                      CPLSPrintf("GDB:PREVIEW:%d:%s",
                                                 iPart+1, osFilename.c_str()));
                aosSubDS.SetNameValue("SUBDATASET_1_DESC", "Preview image");
                poDS->SetMetadata(aosSubDS.List(), "SUBDATASETS");
            }
        }
        else
        {
            CPLStringList aosSubDS;
            for( int i = 0; i < poDS->m_pMPIF->parts(); i++ )
            {
                const auto& header = poDS->m_pMPIF->header(i);
                aosSubDS.SetNameValue(CPLSPrintf("SUBDATASET_%d_NAME", i+1),
                                      CPLSPrintf("GDB:%d:%s", i+1,
                                                 poOpenInfo->pszFilename));
                aosSubDS.SetNameValue(CPLSPrintf("SUBDATASET_%d_DESC", i+1),
                                      header.name().c_str());
            }
            poDS->SetMetadata(aosSubDS.List(), "SUBDATASETS");
        }

        poDS->SetPamFlags(0);

        // Initialize any PAM information.
        poDS->SetDescription(poOpenInfo->pszFilename);
        poDS->TryLoadXML();

        return poDS.release();
    }
    catch( const std::exception& e )
    {
        CPLError(CE_Failure, CPLE_AppDefined, "OpenGDB: %s", e.what());
        return nullptr;
    }
}

/************************************************************************/
/*                          getPixelType()                              */
/************************************************************************/

static PixelType getPixelType(GDALDataType eSrcDT, char ** papszOptions)
{
    PixelType pixelType =
        (eSrcDT == GDT_Byte) ? HALF :
        (eSrcDT == GDT_Int16 ||
            eSrcDT == GDT_UInt16 ||
            eSrcDT == GDT_UInt32) ? UINT: FLOAT;
    const char* pszPixelType =
        CSLFetchNameValueDef(papszOptions, "PIXEL_TYPE", "");
    if( EQUAL(pszPixelType, "HALF") )
        pixelType = HALF;
    else if( EQUAL(pszPixelType, "FLOAT") )
        pixelType = FLOAT;
    else if( EQUAL(pszPixelType, "UINT") )
        pixelType = UINT;
    return pixelType;
}

static void WriteSRSInHeader(Header& header, const OGRSpatialReference* poSRS)
{
    char* pszWKT = nullptr;
    const char* apszOptions[] = { "FORMAT=WKT2_2018", nullptr };
    poSRS->exportToWkt(&pszWKT, apszOptions);
    if( pszWKT )
    {
        header.insert ("gdal:crsWkt", StringAttribute (pszWKT));
        CPLFree(pszWKT);
    }
}

static void WriteGeoTransformInHeader(Header& header, const double* padfGT)
{
    M33d gt;
    gt[0][0] = padfGT[1];
    gt[0][1] = padfGT[2];
    gt[0][2] = padfGT[0];
    gt[1][0] = padfGT[4];
    gt[1][1] = padfGT[5];
    gt[1][2] = padfGT[3];
    gt[2][0] = 0;
    gt[2][1] = 0;
    gt[2][2] = 1;
    header.insert ("gdal:geoTransform", M33dAttribute (gt));
}

static void WriteMetadataInHeader(Header& header, CSLConstList papszMD)
{
    for( CSLConstList papszIter = papszMD; papszIter && *papszIter; ++papszIter )
    {
        char* pszKey = nullptr;
        const char* pszValue = CPLParseNameValue(*papszIter, &pszKey);
        if( pszKey && pszValue )
        {
            header.insert( (std::string("gdal:") + pszKey).c_str(),
                            StringAttribute(pszValue) );
        }
        CPLFree(pszKey);
    }
}

static void FillHeaderFromDataset(Header& header, GDALDataset* poDS)
{
    const auto poSRS = poDS->GetSpatialRef();
    if( poSRS )
    {
        WriteSRSInHeader(header, poSRS);
    }

    double adfGT[6];
    if( poDS->GetGeoTransform(adfGT) == CE_None )
    {
        WriteGeoTransformInHeader(header, adfGT);
    }

    WriteMetadataInHeader(header, poDS->GetMetadata());
}

static void FillHeaderFromOptions(Header& header, CSLConstList papszOptions)
{
    const char* pszDWACompressLevel = CSLFetchNameValue(
        papszOptions, "DWA_COMPRESSION_LEVEL");
    if( pszDWACompressLevel )
    {
        header.insert( "dwaCompressionLevel", FloatAttribute(
            static_cast<float>(CPLAtof(pszDWACompressLevel)) ) );
    }
}

/************************************************************************/
/*                             CreateCopy()                             */
/************************************************************************/

GDALDataset *GDALGDBDataset::CreateCopy( const char* pszFilename,
                                         GDALDataset *poSrcDS,
                                         int, char ** papszOptions,
                                         GDALProgressFunc pfnProgress,
                                         void * pProgressData )
{
    const int nBands = poSrcDS->GetRasterCount();
    const int nXSize = poSrcDS->GetRasterXSize();
    const int nYSize = poSrcDS->GetRasterYSize();
    if( nBands == 0 )
        return nullptr;

    bool bRGB_or_RGBA = false;
    if( (nBands == 3 || nBands == 4) )
    {
        bRGB_or_RGBA = true;
        for( int iBand = 0; iBand < nBands; iBand++ )
        {
            bRGB_or_RGBA &= (poSrcDS->GetRasterBand(iBand+1)->
                        GetColorInterpretation() == GCI_RedBand + iBand);
        }
    }

    const bool bPreview =
        CPLTestBool(CSLFetchNameValueDef(papszOptions, "PREVIEW", "NO")) &&
        (nXSize > 100 || nYSize > 100);
    const GDALDataType eSrcDT = poSrcDS->GetRasterBand(1)->GetRasterDataType();
    if( bPreview && !(bRGB_or_RGBA && eSrcDT == GDT_Byte) )
    {
        CPLError(CE_Failure, CPLE_NotSupported,
                 "Preview creation only supported on RGB/RGBA images of type Byte");
        return nullptr;
    }
    const PixelType pixelType = getPixelType(eSrcDT, papszOptions);
    const bool bRescaleDiv255 =
        pixelType == HALF && bRGB_or_RGBA && eSrcDT == GDT_Byte &&
        CPLTestBool(CSLFetchNameValueDef(papszOptions, "AUTO_RESCALE", "YES"));

    setNumThreads();

    CPLString osTmpOvrFile;
    try
    {
        VSILFILE* fp = VSIFOpenL(pszFilename, "wb+");
        if( fp == nullptr )
        {
            CPLError(CE_Failure, CPLE_FileIO, "Cannot create %s", pszFilename);
            return nullptr;
        }
        GDALGDBIOStream ostream(fp, pszFilename);

        std::vector<std::string> channelNames;
        if( bRGB_or_RGBA )
        {
            channelNames.push_back("R");
            channelNames.push_back("G");
            channelNames.push_back("B");
            if( nBands == 4 )
            {
                channelNames.push_back("A");
            }
        }
        else
        {
            for( int iBand = 0; iBand < nBands; iBand++ )
            {
                channelNames.push_back(CPLSPrintf("Band%d", iBand+1));
            }
        }

        Header header(nXSize, nYSize);

        if( bPreview )
        {
            const int previewWidth = 100;
            const int previewHeight = std::max(1,
                static_cast<int>(static_cast<GIntBig>(previewWidth) * nYSize / nXSize));
            std::vector<PreviewRgba> pixels(previewWidth * previewHeight);
            if( poSrcDS->RasterIO(
                    GF_Read, 0, 0, nXSize, nYSize,
                    &pixels[0],
                    previewWidth, previewHeight,
                    GDT_Byte,
                    nBands, nullptr,
                    4, 4 * previewWidth, 1, nullptr) == CE_None )
            {
                header.setPreviewImage
                    (PreviewImage (previewWidth, previewHeight, &pixels[0]));
            }
        }

        FillHeaderFromDataset(header, poSrcDS);

        const char* pszCompress = CSLFetchNameValueDef(papszOptions, "COMPRESS", "");
        if( pszCompress[0] != '\0' )
        {
            bool bFound = false;
            for( size_t i = 0; i < CPL_ARRAYSIZE(apszCompressions); i++ )
            {
                if( EQUAL(pszCompress, apszCompressions[i]) )
                {
                    bFound = true;
                    header.compression() = static_cast<Compression>(i);
                    break;
                }
            }
            if( !bFound )
            {
                CPLError(CE_Failure, CPLE_AppDefined,
                         "Unknown compression %s", pszCompress);
                return nullptr;
            }
        }

        FillHeaderFromOptions(header, papszOptions);

        std::vector<half> bufferHalf;
        std::vector<float> bufferFloat;
        std::vector<GUInt32> bufferUInt;
        const size_t pixelTypeSize = (pixelType == HALF) ? 2 : 4;
        const GDALDataType eDT = (pixelType == UINT) ? GDT_UInt32 : GDT_Float32;
        const GSpacing nDTSize = GDALGetDataTypeSizeBytes(eDT);

        const bool bTiled = CPLTestBool(CSLFetchNameValueDef(
                                            papszOptions, "TILED", "YES"));

        int nChunkXSize;
        int nChunkYSize;
        const int nBlockXSize = atoi(
            CSLFetchNameValueDef(papszOptions, "BLOCKXSIZE", "256"));
        const int nBlockYSize = atoi(
            CSLFetchNameValueDef(papszOptions, "BLOCKYSIZE", "256"));
        if( nBlockXSize <= 8 || nBlockYSize <= 8 || nBlockXSize >= 8192 ||
            nBlockYSize >= 8192 )
        {
            CPLError(CE_Failure, CPLE_NotSupported,
                     "Invalid block size");
            return nullptr;
        }
        constexpr int MAX_BUFFER_SIZE = 10 * 1024 * 1024;

        const bool bBuildOvr = CPLTestBool(CSLFetchNameValueDef(
                                            papszOptions, "OVERVIEWS", "NO"));
        if( bBuildOvr && !bTiled )
        {
            CPLError(CE_Failure, CPLE_NotSupported,
                     "Overviews only supported on tiled images");
            return nullptr;
        }

        if( bTiled )
        {
            header.setType(TILEDIMAGE);
            header.setTileDescription(TileDescription(nBlockXSize, nBlockYSize,
                                                      bBuildOvr ? MIPMAP_LEVELS : ONE_LEVEL,
                                                      ROUND_UP));
            nChunkYSize = nBlockYSize;
            nChunkXSize = std::min(std::max(nBlockXSize,
                static_cast<int>(MAX_BUFFER_SIZE /
                    (pixelTypeSize * nBands * nBlockYSize) / nBlockXSize * nBlockXSize)),
                nXSize);
        }
        else
        {
            header.setType(SCANLINEIMAGE);
            nChunkXSize = nXSize;
            nChunkYSize = std::min(std::max(1,
                static_cast<int>(MAX_BUFFER_SIZE / (pixelTypeSize * nBands * nXSize))),
                nYSize);
        }
        char* sliceBuffer;
        if( pixelType == UINT )
        {
            bufferUInt.resize(nBands * nChunkXSize * nChunkYSize);
            sliceBuffer = reinterpret_cast<char*>(bufferUInt.data());
        }
        else
        {
            bufferFloat.resize(nBands * nChunkXSize * nChunkYSize);
            if( pixelType == HALF )
            {
                bufferHalf.resize(nBands * nChunkXSize * nChunkYSize);
                sliceBuffer = reinterpret_cast<char*>(bufferHalf.data());
            }
            else
            {
                sliceBuffer = reinterpret_cast<char*>(bufferFloat.data());
            }
        }

        for( const auto& channelName: channelNames )
        {
            header.channels().insert(channelName, Channel(pixelType));
        }

        MultiPartOutputFile mpof (ostream, &header, 1);
        if( bTiled )
        {
            TiledOutputPart op(mpof, 0);

            if( bBuildOvr )
            {
                if( nBlockXSize != nBlockYSize )
                {
                    CPLError(CE_Failure, CPLE_NotSupported,
                             "Overview building only works if BLOCKXSIZE=BLOCKYSIZE");
                    return nullptr;
                }
                if( nBlockXSize < 64 || nBlockXSize > 4096 ||
                    !CPLIsPowerOfTwo(nBlockXSize) )
                {
                    CPLError(CE_Failure, CPLE_NotSupported,
                             "Overview building only works if "
                             "BLOCKXSIZE=BLOCKYSIZE is a power of 2 "
                             "between 64 and 4096.");
                    return nullptr;
                }
            }

            const auto writeTiles = [nChunkXSize, nChunkYSize,
                                     nBlockXSize, nBlockYSize,
                                     nBands,
                                     pixelType,
                                     pixelTypeSize,
                                     sliceBuffer,
                                     eDT,
                                     nDTSize,
                                     bRescaleDiv255,
                                     &channelNames,
                                     &op,
                                     &bufferFloat,
                                     &bufferHalf,
                                     &bufferUInt]
                (GDALDataset* l_poDS,
                 int iLevel,
                 GDALProgressFunc l_pfnProgress,
                 void * l_pProgressData)
            {
                const int l_nXSize = l_poDS->GetRasterXSize();
                const int l_nYSize = l_poDS->GetRasterYSize();
                const int nXBlocks = DIV_ROUND_UP(l_nXSize, nBlockXSize);
                const int nYBlocks = DIV_ROUND_UP(l_nYSize, nBlockYSize);
                for( int y = 0; y < l_nYSize; y += nChunkYSize)
                {
                    const int nLinesToRead = std::min(nChunkYSize, l_nYSize - y);
                    for( int x = 0; x < l_nXSize; x += nChunkXSize)
                    {
                        const int nColsToRead = std::min(nChunkXSize, l_nXSize - x);
                        FrameBuffer fb;
                        for( int iBand = 0; iBand < nBands; iBand++ )
                        {
                            const auto slice =
                                Slice(pixelType,
                                    sliceBuffer +
                                        iBand * pixelTypeSize * nChunkXSize * nChunkYSize -
                                        (x * pixelTypeSize + y * pixelTypeSize * nChunkXSize),
                                    pixelTypeSize, pixelTypeSize * nChunkXSize);
                            fb.insert(channelNames[iBand], slice);
                        }
                        if( l_poDS->RasterIO(
                            GF_Read,
                            x, y, nColsToRead, nLinesToRead,
                            !bufferFloat.empty() ?
                                reinterpret_cast<GByte*>(&bufferFloat[0]):
                                reinterpret_cast<GByte*>(&bufferUInt[0]),
                            nColsToRead, nLinesToRead, eDT,
                            nBands, nullptr,
                            nDTSize,
                            nDTSize * nChunkXSize,
                            nDTSize * nChunkXSize * nChunkYSize,
                            nullptr) != CE_None )
                        {
                            return false;
                        }
                        if( pixelType == HALF )
                        {
                            const size_t nPixelsInBuffer =
                                static_cast<size_t>(nChunkXSize) * nChunkYSize * nBands;
                            if( bRescaleDiv255 )
                            {
                                for( size_t i = 0; i < nPixelsInBuffer; i++ )
                                {
                                    bufferHalf[i] = bufferFloat[i] / 255.0f;
                                }
                            }
                            else
                            {
                                for( size_t i = 0; i < nPixelsInBuffer; i++ )
                                {
                                    bufferHalf[i] = bufferFloat[i];
                                }
                            }
                        }
                        op.setFrameBuffer(fb);
                        const int blockXMax = (x + nColsToRead - 1) / nBlockXSize;
                        const int blockYMax = (y + nLinesToRead - 1) / nBlockYSize;
                        op.writeTiles(x / nBlockXSize, blockXMax,
                                      y / nBlockYSize, blockYMax,
                                      iLevel);
                        if( l_pfnProgress &&
                            !l_pfnProgress(
                                (static_cast<double>(blockYMax) * nXBlocks +
                                blockXMax + 1) / nXBlocks / nYBlocks , "",
                                        l_pProgressData) )
                        {
                            return false;
                        }
                    }
                }
                return true;
            };

            struct ScaledProgressReleaser
            {
                void operator()(void* progress) const {
                    GDALDestroyScaledProgress(progress); }
            };

            using ScaledProgressUniquePtr = std::unique_ptr<void, ScaledProgressReleaser>;
            ScaledProgressUniquePtr progress;

            // Write full resolution imagery
            if( bBuildOvr )
                progress.reset(GDALCreateScaledProgress(0, 0.5, pfnProgress, pProgressData));
            else
                progress.reset(GDALCreateScaledProgress(0, 1, pfnProgress, pProgressData));
            if( !writeTiles(poSrcDS, 0, GDALScaledProgress, progress.get()) )
            {
                if( !osTmpOvrFile.empty() )
                    VSIUnlink(osTmpOvrFile);
                return nullptr;
            }

            if( bBuildOvr )
            {
                // First build overviews in a temporary GTiff file
                GDALDefaultOverviews oOvr;
                oOvr.Initialize(poSrcDS);
                std::vector<int> anOvrFactors;
                for( int i = 1; i < op.numLevels(); i++ )
                    anOvrFactors.push_back(1 << i);
                std::vector<int> anBands;
                for( int iBand = 0; iBand < nBands; iBand++ )
                    anBands.push_back(iBand+1);
                CPLSetThreadLocalConfigOption("GDAL_TIFF_OVR_BLOCKSIZE",
                                              CPLSPrintf("%d", nBlockXSize));
                const CPLString osTmpOvrFileRadix(CPLSPrintf("%s_tmp",pszFilename));
                osTmpOvrFile = osTmpOvrFileRadix + ".ovr";
                progress.reset(GDALCreateScaledProgress(0.5, 0.8, pfnProgress, pProgressData));
                if( oOvr.BuildOverviews(osTmpOvrFileRadix,
                                    CSLFetchNameValueDef(papszOptions,
                                        "OVERVIEW_RESAMPLING", "CUBIC"),
                                    static_cast<int>(anOvrFactors.size()),
                                    &anOvrFactors[0],
                                    nBands, &anBands[0],
                                    GDALScaledProgress, progress.get()) != CE_None )
                {
                    CPLSetThreadLocalConfigOption("GDAL_TIFF_OVR_BLOCKSIZE",
                                                  nullptr);
                    VSIUnlink(osTmpOvrFile);
                    return nullptr;
                }
                CPLSetThreadLocalConfigOption("GDAL_TIFF_OVR_BLOCKSIZE",
                                              nullptr);

                // Transfer overviews from temporary file to main image
                std::unique_ptr<GDALDataset> poOvrDS(GDALDataset::Open(osTmpOvrFile));
                if( !poOvrDS )
                    return nullptr;
                const int nOvrs = 1 + poOvrDS->GetRasterBand(1)->GetOverviewCount();
                for( int i = 0; i < nOvrs; i++ )
                {
                    auto poThisOvrDS = (i == 0) ? poOvrDS.get() :
                        poOvrDS->GetRasterBand(1)->GetOverview(i-1)->GetDataset();
                    CPLAssert(poThisOvrDS);
                    if( i == 0 )
                        progress.reset(GDALCreateScaledProgress(
                                       0.8, nOvrs == 1 ? 1.0 : 0.9,
                                       pfnProgress, pProgressData));
                    else if( i == 1 )
                        progress.reset(GDALCreateScaledProgress(
                                       0.9, nOvrs == 2 ? 1.0 : 0.95,
                                       pfnProgress, pProgressData));
                    else
                        progress.reset(GDALCreateScaledProgress(
                                       0.95 + 0.05 * (i - 2) / (nOvrs - 2),
                                       0.95 + 0.05 * (i - 2 + 1) / (nOvrs - 2),
                                       pfnProgress, pProgressData));
                    if( !writeTiles(poThisOvrDS, i+1, GDALScaledProgress,
                                    progress.get()) )
                    {
                        poOvrDS.reset();
                        VSIUnlink(osTmpOvrFile);
                        return nullptr;
                    }
                }

                poOvrDS.reset();
                VSIUnlink(osTmpOvrFile);
            }
        }
        else
        {
            OutputPart op(mpof, 0);

            for( int y = 0; y < nYSize; y+= nChunkYSize)
            {
                FrameBuffer fb;
                const int nLinesToRead = std::min(nChunkYSize, nYSize - y);
                for( int iBand = 0; iBand < nBands; iBand++ )
                {
                    const auto slice =
                        Slice(pixelType,
                            sliceBuffer +
                                iBand * pixelTypeSize * nXSize * nLinesToRead -
                                y * pixelTypeSize * nXSize,
                            pixelTypeSize, pixelTypeSize * nXSize);
                    fb.insert(channelNames[iBand], slice);
                }
                if( poSrcDS->RasterIO(
                    GF_Read,
                    0, y, nXSize, nLinesToRead,
                    !bufferFloat.empty() ?
                        reinterpret_cast<GByte*>(&bufferFloat[0]):
                        reinterpret_cast<GByte*>(&bufferUInt[0]),
                    nXSize, nLinesToRead, eDT,
                    nBands, nullptr,
                    nDTSize,
                    nDTSize * nXSize,
                    nDTSize * nXSize * nLinesToRead,
                    nullptr) != CE_None )
                {
                    return nullptr;
                }
                if( pixelType == HALF )
                {
                    for( size_t i = 0;
                            i < static_cast<size_t>(nXSize) * nLinesToRead * nBands; i++ )
                    {
                        // cppcheck-suppress unreadVariable
                        bufferHalf[i] = bufferFloat[i];
                    }
                }
                op.setFrameBuffer(fb);
                op.writePixels(nLinesToRead);
                if( pfnProgress &&
                    !pfnProgress(static_cast<double>(y+nLinesToRead) / nYSize, "",
                                pProgressData) )
                {
                    return nullptr;
                }
            }
        }
    }
    catch( const std::exception& e )
    {
        if( !osTmpOvrFile.empty() )
            VSIUnlink(osTmpOvrFile);
        CPLError(CE_Failure, CPLE_AppDefined, "OpenGDB: %s", e.what());
        return nullptr;
    }
    GDALOpenInfo oOpenInfo(pszFilename, GA_ReadOnly);
    return GDALGDBDataset::Open(&oOpenInfo);
}

/************************************************************************/
/*                         GDALGDBWritableDataset                       */
/************************************************************************/

class GDALGDBWritableDataset final: public GDALPamDataset
{
    friend class GDALGDBDataset;
    friend class GDALGDBWritableRasterBand;

    PixelType m_pixelType = HALF;
    int m_nBlockXSize = 0;
    int m_nBlockYSize = 0;

    // Keep stream before others, so that it is destroyed last
    std::unique_ptr<OStream> m_pOStream{};

    std::unique_ptr<TiledOutputPart> m_pTOP{};
    std::unique_ptr<MultiPartOutputFile> m_pMPOF{};

    std::vector<std::string> m_channelNames{};

    bool m_bTriedWritingHeader = false;
    std::vector<half> m_bufferHalf{};
    std::vector<float> m_bufferFloat{};
    std::vector<GUInt32> m_bufferUInt{};
    size_t m_nBufferEltSize = 0;
    char* m_pSliceBuffer = nullptr;

    OGRSpatialReference m_oSRS{};
    double m_adfGT[6] = {0,1,0,0,0,1};
    bool m_bHasGT = false;

    CPLStringList m_aosMetadata{};

    std::vector<bool> m_abWrittenBlocks{};
    size_t m_nXBlocks = 0;

    bool m_bRescaleDiv255 = false;

    Header m_header;

    void WriteHeader();

public:
        GDALGDBWritableDataset(int nXSize, int nYSize): m_header(nXSize, nYSize)
        {
            nRasterXSize = nXSize;
            nRasterYSize = nYSize;
        }
        ~GDALGDBWritableDataset() override;

        CPLErr SetGeoTransform(double* adfGT) override;
        CPLErr SetSpatialRef(const OGRSpatialReference* poSRS) override;

        const OGRSpatialReference* GetSpatialRef() const override;
        CPLErr GetGeoTransform(double* adfGT) override;

        CPLErr  SetMetadata( char **, const char * = "" ) override;
        CPLErr  SetMetadataItem( const char*, const char*,
                                     const char* = "" ) override;

        char** GetMetadata( const char* pszDomain = "" ) override;
        const char* GetMetadataItem( const char* pszName,
                                     const char* pszDomain = "" ) override;
};

/************************************************************************/
/*                       ~GDALGDBWritableDataset()                      */
/************************************************************************/

GDALGDBWritableDataset::~GDALGDBWritableDataset()
{
    WriteHeader();
    FlushCache();
}

/************************************************************************/
/*                            SetGeoTransform()                         */
/************************************************************************/

CPLErr GDALGDBWritableDataset::SetGeoTransform(double* adfGT)
{
    if( m_bTriedWritingHeader )
    {
        CPLError(CE_Warning, CPLE_AppDefined,
                 "SetGeoTransform() called after writing pixels. Will go to PAM");
        return GDALPamDataset::SetGeoTransform(adfGT);
    }
    m_bHasGT = true;
    memcpy(m_adfGT, adfGT, 6 * sizeof(double));
    return CE_None;
}

/************************************************************************/
/*                            SetSpatialRef()                           */
/************************************************************************/

CPLErr GDALGDBWritableDataset::SetSpatialRef(const OGRSpatialReference* poSRS)
{
    if( m_bTriedWritingHeader )
    {
        CPLError(CE_Warning, CPLE_AppDefined,
                 "SetSpatialRef() called after writing pixels. Will go to PAM");
        return GDALPamDataset::SetSpatialRef(poSRS);
    }
    if( poSRS )
        m_oSRS = *poSRS;
    else
        m_oSRS.Clear();
    return CE_None;
}

/************************************************************************/
/*                             SetMetadata()                            */
/************************************************************************/

CPLErr GDALGDBWritableDataset::SetMetadata( char ** papszMD,
                                            const char* pszDomain)
{
    if( pszDomain == nullptr || pszDomain[0] == 0 )
    {
        m_aosMetadata = CSLDuplicate(papszMD);
        if( m_bTriedWritingHeader )
        {
            CPLError(CE_Warning, CPLE_AppDefined,
                     "SetMetadata() called after writing pixels. Will go to PAM");
        }
        else
        {
            return CE_None;
        }
    }
    return GDALPamDataset::SetMetadata(papszMD, pszDomain);
}

/************************************************************************/
/*                           SetMetadataItem()                          */
/************************************************************************/

CPLErr GDALGDBWritableDataset::SetMetadataItem( const char* pszName,
                                                const char* pszValue,
                                                const char* pszDomain)
{
    if( pszDomain == nullptr || pszDomain[0] == 0 )
    {
        m_aosMetadata.SetNameValue(pszName, pszValue);
        if( m_bTriedWritingHeader )
        {
            CPLError(CE_Warning, CPLE_AppDefined,
                     "SetMetadata() called after writing pixels. Will go to PAM");
        }
        else
        {
            return CE_None;
        }
    }
    return GDALPamDataset::SetMetadataItem(pszName, pszValue, pszDomain);
}

/************************************************************************/
/*                             GetMetadata()                            */
/************************************************************************/

char** GDALGDBWritableDataset::GetMetadata( const char* pszDomain )
{
    if( pszDomain == nullptr || pszDomain[0] == 0 )
    {
        return m_aosMetadata.List();
    }
    return GDALPamDataset::GetMetadata(pszDomain);
}

/************************************************************************/
/*                           GetMetadataItem()                          */
/************************************************************************/

const char* GDALGDBWritableDataset::GetMetadataItem( const char* pszName,
                                                     const char* pszDomain )
{
    if( pszDomain == nullptr || pszDomain[0] == 0 )
    {
        return m_aosMetadata.FetchNameValue(pszName);
    }
    return GDALPamDataset::GetMetadataItem(pszName, pszDomain);
}

/************************************************************************/
/*                          GetSpatialRef()                             */
/************************************************************************/

const OGRSpatialReference* GDALGDBWritableDataset::GetSpatialRef() const
{
    const auto* poPamSRS = GDALPamDataset::GetSpatialRef();
    if( poPamSRS )
        return poPamSRS;
    return m_oSRS.IsEmpty() ? nullptr : &m_oSRS;
}

/************************************************************************/
/*                         GetGeoTransform()                            */
/************************************************************************/

CPLErr GDALGDBWritableDataset::GetGeoTransform(double* adfGT)
{
    if( GDALPamDataset::GetGeoTransform(adfGT) == CE_None )
    {
        return CE_None;
    }
    memcpy(adfGT, m_adfGT, 6 * sizeof(double));
    return m_bHasGT ? CE_None : CE_Failure;
}

/************************************************************************/
/*                             WriteHeader()                            */
/************************************************************************/

void GDALGDBWritableDataset::WriteHeader()
{
    if( m_bTriedWritingHeader )
        return;
    m_bTriedWritingHeader = true;

    try
    {
        FillHeaderFromDataset(m_header, this);

        bool bRGB_or_RGBA = false;
        if( nBands == 3 || nBands == 4 )
        {
            bRGB_or_RGBA = true;
            for( int i = 0; i < nBands; i++ )
            {
                bRGB_or_RGBA &=
                    GetRasterBand(i+1)->GetColorInterpretation() == GCI_RedBand + i;
            }
        }
        m_bRescaleDiv255 &=
            m_pixelType == HALF && bRGB_or_RGBA &&
            GetRasterBand(1)->GetRasterDataType() == GDT_Byte;

        if( bRGB_or_RGBA )
        {
            m_channelNames.push_back("R");
            m_channelNames.push_back("G");
            m_channelNames.push_back("B");
            if( nBands == 4 )
            {
                m_channelNames.push_back("A");
            }
        }
        else
        {
            for( int iBand = 0; iBand < nBands; iBand++ )
            {
                m_channelNames.push_back(CPLSPrintf("Band%d", iBand+1));
            }
        }

        for( int i = 0; i < nBands; i++ )
        {
            m_header.channels().insert(m_channelNames[i], Channel(m_pixelType));
        }

        m_pMPOF.reset(new MultiPartOutputFile(*m_pOStream, &m_header, 1));
        m_pTOP.reset(new TiledOutputPart(*m_pMPOF, 0));

        const size_t nElts =
            static_cast<size_t>(nBands) * m_nBlockXSize * m_nBlockYSize;
        if( m_pixelType == HALF )
        {
            m_bufferHalf.resize(nElts);
            m_bufferFloat.resize(nElts / nBands);
            m_pSliceBuffer = reinterpret_cast<char*>(&m_bufferHalf[0]);
            m_nBufferEltSize = sizeof(half);
        }
        else if( m_pixelType == FLOAT )
        {
            m_bufferFloat.resize(nElts);
            m_pSliceBuffer = reinterpret_cast<char*>(&m_bufferFloat[0]);
            m_nBufferEltSize = sizeof(float);
        }
        else
        {
            m_bufferUInt.resize(nElts);
            m_pSliceBuffer = reinterpret_cast<char*>(&m_bufferUInt[0]);
            m_nBufferEltSize = sizeof(unsigned int);
        }
    }
    catch( const std::exception& e )
    {
        CPLError(CE_Failure, CPLE_AppDefined, "OpenGDB: %s", e.what());
        m_pTOP.reset();
        m_pMPOF.reset();
    }
}

/************************************************************************/
/*                       GDALGDBWritableRasterBand                      */
/************************************************************************/

class GDALGDBWritableRasterBand final: public GDALPamRasterBand
{
        GDALColorInterp m_eInterp = GCI_Undefined;

    protected:
        CPLErr IReadBlock(int, int, void*) override;
        CPLErr IWriteBlock(int, int, void*) override;

    public:
        GDALGDBWritableRasterBand(GDALGDBWritableDataset* poDSIn,
                                  int nBandIn,
                                  GDALDataType eTypeIn);

        CPLErr SetColorInterpretation(GDALColorInterp eInterp) override
            { m_eInterp = eInterp; return CE_None; }
        GDALColorInterp GetColorInterpretation() override
            { return m_eInterp; }
};

/************************************************************************/
/*                       GDALGDBWritableRasterBand()                    */
/************************************************************************/

GDALGDBWritableRasterBand::GDALGDBWritableRasterBand(
                                    GDALGDBWritableDataset* poDSIn,
                                    int nBandIn,
                                    GDALDataType eTypeIn)
{
    poDS = poDSIn;
    nBand = nBandIn;
    nRasterXSize = poDSIn->GetRasterXSize();
    nRasterYSize = poDSIn->GetRasterYSize();
    nBlockXSize = poDSIn->m_nBlockXSize;
    nBlockYSize = poDSIn->m_nBlockYSize;
    eDataType = eTypeIn;
}

/************************************************************************/
/*                           IReadBlock()                               */
/************************************************************************/

CPLErr GDALGDBWritableRasterBand::IReadBlock(int nBlockXOff,
                                             int nBlockYOff,
                                             void* pImage)
{
    auto poGDS = cpl::down_cast<GDALGDBWritableDataset*>(poDS);
    if( !poGDS->m_abWrittenBlocks[nBlockYOff * poGDS->m_nXBlocks + nBlockXOff] )
    {
        const size_t nPixelsInBlock =
            static_cast<size_t>(nBlockXSize) * nBlockYSize;
        memset(pImage, 0, nPixelsInBlock * GDALGetDataTypeSizeBytes(eDataType));
        return CE_None;
    }
    CPLError(CE_Failure, CPLE_AppDefined,
             "Reading blocks in a GDB dataset created by Create() is not "
             "supported");
    return CE_Failure;
}

/************************************************************************/
/*                           IWriteBlock()                              */
/************************************************************************/

CPLErr GDALGDBWritableRasterBand::IWriteBlock(int nBlockXOff,
                                              int nBlockYOff,
                                              void* pImage)
{
    auto poGDS = cpl::down_cast<GDALGDBWritableDataset*>(poDS);
    poGDS->WriteHeader();
    if( !poGDS->m_pTOP )
        return CE_Failure;

    poGDS->m_abWrittenBlocks[nBlockYOff * poGDS->m_nXBlocks + nBlockXOff] = true;

    bool bAllBlocksDirty = true;
    std::vector<GDALRasterBlock*> apoBlocks;
    apoBlocks.resize(poGDS->nBands);
    for( int iBand = 0; iBand < poGDS->nBands; ++iBand )
    {
        if( iBand + 1 != nBand )
        {
            apoBlocks[iBand] =
                cpl::down_cast<GDALGDBWritableRasterBand *>(
                    poGDS->GetRasterBand( iBand + 1 ))
                        ->TryGetLockedBlockRef( nBlockXOff, nBlockYOff );

            if( apoBlocks[iBand] == nullptr )
            {
                bAllBlocksDirty = false;
                break;
            }
            else if( !apoBlocks[iBand]->GetDirty() )
            {
                apoBlocks[iBand]->DropLock();
                apoBlocks[iBand] = nullptr;
                bAllBlocksDirty = false;
                break;
            }
        }
        else
        {
            apoBlocks[iBand] = nullptr;
        }
    }
    if( !bAllBlocksDirty )
    {
        CPLError(CE_Warning, CPLE_AppDefined,
                 "For block (%d, %d), blocks for some bands are not available "
                 "in the cache. Corresponding data will be assumed to be zero.",
                 nBlockXOff, nBlockYOff);
    }

    CPLErr eErr = CE_None;
    try
    {
        FrameBuffer fb;
        const int x = nBlockXOff * nBlockXSize;
        const int y = nBlockYOff * nBlockYSize;
        const size_t nPixelsInBlock =
            static_cast<size_t>(nBlockXSize) * nBlockYSize;
        const GDALDataType eDstDT =
            poGDS->m_pixelType == UINT ? GDT_UInt32 : GDT_Float32;
        for( int iBand = 0; iBand < poGDS->nBands; iBand++ )
        {
            char* const dstPtr = poGDS->m_pSliceBuffer +
                iBand * poGDS->m_nBufferEltSize * nPixelsInBlock;
            const auto slice =
                Slice(poGDS->m_pixelType,
                      dstPtr -
                        (x * poGDS->m_nBufferEltSize +
                         y * poGDS->m_nBufferEltSize * nBlockXSize),
                      poGDS->m_nBufferEltSize,
                      poGDS->m_nBufferEltSize * nBlockXSize);
            fb.insert(poGDS->m_channelNames[iBand], slice);

            const void* srcPtr = nullptr;
            if( iBand+1 == nBand)
                srcPtr = pImage;
            else if( apoBlocks[iBand] )
                srcPtr = apoBlocks[iBand]->GetDataRef();
            else
            {
                memset(poGDS->m_pSliceBuffer +
                            iBand * poGDS->m_nBufferEltSize * nPixelsInBlock,
                       0,
                       nPixelsInBlock * poGDS->m_nBufferEltSize);
                continue;
            }

            GDALCopyWords64(
                srcPtr,
                eDataType,
                GDALGetDataTypeSizeBytes(eDataType),
                poGDS->m_pixelType == HALF ?
                    static_cast<void*>(&poGDS->m_bufferFloat[0]):
                    static_cast<void*>(dstPtr),
                eDstDT,
                GDALGetDataTypeSizeBytes(eDstDT),
                static_cast<GPtrDiff_t>(nPixelsInBlock));
            if( poGDS->m_pixelType == HALF )
            {
                if( poGDS->m_bRescaleDiv255 )
                {
                    for( size_t i = 0; i < nPixelsInBlock; i++ )
                    {
                        poGDS->m_bufferHalf[iBand * nPixelsInBlock + i] =
                            poGDS->m_bufferFloat[i] / 255.0f;
                    }
                }
                else
                {
                    for( size_t i = 0; i < nPixelsInBlock; i++ )
                    {
                        poGDS->m_bufferHalf[iBand * nPixelsInBlock + i] =
                            poGDS->m_bufferFloat[i];
                    }
                }
            }
        }

        poGDS->m_pTOP->setFrameBuffer(fb);
        poGDS->m_pTOP->writeTile(nBlockXOff, nBlockYOff);
    }
    catch( const std::exception& e )
    {
        CPLError(CE_Failure, CPLE_AppDefined, "OpenGDB: %s", e.what());
        eErr = CE_Failure;
    }

    for( int iBand = 0; iBand < poGDS->nBands; ++iBand )
    {
        if( apoBlocks[iBand] )
        {
            apoBlocks[iBand]->MarkClean();
            apoBlocks[iBand]->DropLock();
        }
    }

    return eErr;
}

/************************************************************************/
/*                            Create()                                  */
/************************************************************************/

GDALDataset *GDALGDBDataset::Create( const char * pszFilename,
                                     int nXSize, int nYSize, int nBands,
                                     GDALDataType eType, char ** papszOptions )
{
    if( nBands == 0 )
        return nullptr;
    const PixelType pixelType = getPixelType(eType, papszOptions);

    if( !CPLTestBool(CSLFetchNameValueDef(papszOptions, "TILED", "YES")) )
    {
        CPLError(CE_Failure, CPLE_NotSupported,
                 "Create() only supports tiled mode");
        return nullptr;
    }

    if( CPLTestBool(CSLFetchNameValueDef(papszOptions, "OVERVIEWS", "NO")) )
    {
        CPLError(CE_Failure, CPLE_NotSupported,
                 "Create() does not support overview creation.");
        return nullptr;
    }

    if( CPLTestBool(CSLFetchNameValueDef(papszOptions, "PREVIEW", "NO")) )
    {
        CPLError(CE_Failure, CPLE_NotSupported,
                 "Create() does not support preview creation.");
        return nullptr;
    }

    Compression compression = ZIP_COMPRESSION;
    const char* pszCompress = CSLFetchNameValueDef(papszOptions, "COMPRESS", "");
    if( pszCompress[0] != '\0' )
    {
        bool bFound = false;
        for( size_t i = 0; i < CPL_ARRAYSIZE(apszCompressions); i++ )
        {
            if( EQUAL(pszCompress, apszCompressions[i]) )
            {
                bFound = true;
                compression = static_cast<Compression>(i);
                break;
            }
        }
        if( !bFound )
        {
            CPLError(CE_Failure, CPLE_AppDefined,
                        "Unknown compression %s", pszCompress);
            return nullptr;
        }
    }

    const int nBlockXSize = atoi(
        CSLFetchNameValueDef(papszOptions, "BLOCKXSIZE", "256"));
    const int nBlockYSize = atoi(
        CSLFetchNameValueDef(papszOptions, "BLOCKYSIZE", "256"));
    if( nBlockXSize <= 8 || nBlockYSize <= 8 || nBlockXSize >= 8192 ||
        nBlockYSize >= 8192 )
    {
        CPLError(CE_Failure, CPLE_NotSupported,
                    "Invalid block size");
        return nullptr;
    }

    VSILFILE* fp = VSIFOpenL(pszFilename, "wb+");
    if( fp == nullptr )
    {
        CPLError(CE_Failure, CPLE_FileIO, "Cannot create %s", pszFilename);
        return nullptr;
    }
    auto poDS = std::unique_ptr<GDALGDBWritableDataset>(
                            new GDALGDBWritableDataset(nXSize, nYSize));
    poDS->m_pOStream.reset(new GDALGDBIOStream(fp, pszFilename));
    poDS->eAccess = GA_Update;
    poDS->m_pixelType = pixelType;
    poDS->m_header.compression() = compression;
    poDS->m_header.setType(TILEDIMAGE);
    poDS->m_header.setTileDescription(TileDescription(nBlockXSize, nBlockYSize));
    FillHeaderFromOptions(poDS->m_header, papszOptions);
    poDS->m_nBlockXSize = nBlockXSize;
    poDS->m_nBlockYSize = nBlockYSize;
    poDS->m_nXBlocks = static_cast<size_t>(DIV_ROUND_UP(nXSize, nBlockXSize));
    const size_t nYBlocks = static_cast<size_t>(DIV_ROUND_UP(nYSize, nBlockYSize));
    if( poDS->m_nXBlocks > std::numeric_limits<size_t>::max() / nYBlocks )
    {
        return nullptr;
    }
    try
    {
        poDS->m_abWrittenBlocks.resize(
            poDS->m_nXBlocks * nYBlocks);
    }
    catch( const std::exception& e )
    {
        CPLError(CE_Failure, CPLE_AppDefined, "%s", e.what());
        return nullptr;
    }
    poDS->m_bRescaleDiv255 =
        CPLTestBool(CSLFetchNameValueDef(papszOptions, "AUTO_RESCALE", "YES"));

    if( nBands > 1 )
    {
        poDS->GDALDataset::SetMetadataItem("INTERLEAVE", "PIXEL",
                                           "IMAGE_STRUCTURE");
    }
    for(int i = 0; i < nBands; i++ )
    {
        poDS->SetBand(i+1, new GDALGDBWritableRasterBand(poDS.get(), i+1, eType));
    }
    poDS->SetDescription(pszFilename);
    poDS->TryLoadXML();
    return poDS.release();
}

/************************************************************************/
/*                          GDALRegister_GDB()                          */
/************************************************************************/

void GDALRegister_GDB()

{
    if( !GDAL_CHECK_VERSION("GDB driver") )
        return;

    if( GDALGetDriverByName("GDB") != nullptr )
        return;

    GDALDriver *poDriver = new GDALDriver();

    poDriver->SetDescription("GDB");
    poDriver->SetMetadataItem(GDAL_DCAP_RASTER, "YES");
    poDriver->SetMetadataItem(GDAL_DMD_LONGNAME,
                              "Extended Dynamic Range Image File Format");
    poDriver->SetMetadataItem(GDAL_DMD_HELPTOPIC, "drivers/raster/GDB.html");
    poDriver->SetMetadataItem(GDAL_DMD_EXTENSION, "GDB");
    poDriver->SetMetadataItem(GDAL_DMD_CREATIONOPTIONLIST,
"<CreationOptionList>"
"   <Option name='COMPRESS' type='string-select' default='ZIP'>"
"     <Value>NONE</Value>"
"     <Value>RLE</Value>"
"     <Value>ZIPS</Value>"
"     <Value>ZIP</Value>"
"     <Value>PIZ</Value>"
"     <Value>PXR24</Value>"
"     <Value>B44</Value>"
"     <Value>B44A</Value>"
"     <Value>DWAA</Value>"
"     <Value>DWAB</Value>"
"   </Option>"
"   <Option name='PIXEL_TYPE' type='string-select'>"
"     <Value>HALF</Value>"
"     <Value>FLOAT</Value>"
"     <Value>UINT</Value>"
"   </Option>"
"   <Option name='TILED' type='boolean' description='Use tiling' default='YES'/>"
"   <Option name='BLOCKXSIZE' type='int' description='Tile width' default='256'/>"
"   <Option name='BLOCKYSIZE' type='int' description='Tile height' default='256'/>"
"   <Option name='OVERVIEWS' type='boolean' description='Whether to create overviews' default='NO'/>"
"   <Option name='OVERVIEW_RESAMPLING' type='string' description='Resampling method' default='CUBIC'/>"
"   <Option name='PREVIEW' type='boolean' description='Create a preview' default='NO'/>"
"   <Option name='AUTO_RESCALE' type='boolean' description='Whether to rescale Byte RGB(A) values to 0-1' default='YES'/>"
"   <Option name='DWA_COMPRESSION_LEVEL' type='int' description='DWA compression level'/>"
"</CreationOptionList>"
                              );
    poDriver->SetMetadataItem(GDAL_DMD_SUBDATASETS, "YES");
    poDriver->SetMetadataItem( GDAL_DCAP_VIRTUALIO, "YES" );

    poDriver->pfnOpen = GDALGDBDataset::Open;
    poDriver->pfnIdentify = GDALGDBDataset::Identify;
    poDriver->pfnCreateCopy = GDALGDBDataset::CreateCopy;
    poDriver->pfnCreate = GDALGDBDataset::Create;

    GetGDALDriverManager()->RegisterDriver(poDriver);
}