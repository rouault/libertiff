// SPDX-License-Identifier: MIT
// Copyright 2024, Even Rouault <even.rouault at spatialys.com>

// Canonical URL: https://github.com/libertiff/libertiff/blob/master/libertiff.hpp

#ifndef LIBERTIFF_HPP_INCLUDED
#define LIBERTIFF_HPP_INCLUDED

//////////////////////////////////////////////////////////////
// libertiff = libre TIFF or LIB E(ven) R(ouault) TIFF... ? //
//////////////////////////////////////////////////////////////

#if __cplusplus >= 202002L
#include <bit>  // std::endian
#endif

#include <array>
#include <cassert>
#include <cstring>
#include <limits>
#include <memory>
#include <set>
#include <string>
#include <type_traits>
#include <vector>

#ifndef LIBERTIFF_NS
#define LIBERTIFF_NS libertiff
#endif

/** Libertiff is a C++11 simple, header-only, TIFF reader. It is MIT licensed.
 *
 * Handles both ClassicTIFF and BigTIFF, little-endian or big-endian ordered.
 * 
 * The library does not (yet?) offer codec facilities. It is mostly aimed at
 * browsing through the linked chain of Image File Descriptors and their tags.
 *
 * The library is thread-safe (that is the instances that it returns can
 * be used from multiple threads), if passed FileReader instances are themselves
 * thread-safe.
 *
 * The library does not throw exceptions (but underlying std library might
 * throw exceptions in case of out-of-memory when allocating the tag array)
 *
 * Optional features:
 * - define LIBERTIFF_C_FILE_READER before including that file, so that
 *   the libertiff::CFileReader class is available
 */
namespace LIBERTIFF_NS
{

#if __cplusplus >= 201703L
#define LIBERTIFF_STATIC_ASSERT(x) static_assert(x)
#define LIBERTIFF_CONSTEXPR constexpr
#else
#define LIBERTIFF_STATIC_ASSERT(x) static_assert((x), #x)
#define LIBERTIFF_CONSTEXPR
#endif

template <typename T, typename... Args>
std::unique_ptr<T> make_unique(Args &&...args)
{
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

/** Returns whether the host is little-endian ordered */
inline bool isHostLittleEndian()
{
#if __cplusplus >= 202002L
    return std::endian::native == std::endian::little;
#elif (defined(__BYTE_ORDER__) && defined(__ORDER_LITTLE_ENDIAN__) &&          \
       (__BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__)) ||                         \
    defined(_MSC_VER)
    return true;
#elif defined(__BYTE_ORDER__) && defined(__ORDER_BIG_ENDIAN__) &&              \
    (__BYTE_ORDER__ == __ORDER_BIG_ENDIAN__)
    return false;
#else
    uint32_t one = 1;
    char one_as_char_array[sizeof(uint32_t)];
    std::memcpy(one_as_char_array, &one, sizeof(uint32_t));
    return one_as_char_array[0] == 1;
#endif
}

/** Byte-swap */
template <class T> inline T byteSwap(T v);

/** Byte-swap a uint8_t */
template <> uint8_t byteSwap(uint8_t v)
{
    return v;
}

/** Byte-swap a int8_t */
template <> int8_t byteSwap(int8_t v)
{
    return v;
}

/** Byte-swap a uint16_t */
template <> uint16_t byteSwap(uint16_t v)
{
    return uint16_t((v >> 8) | ((v & 0xff) << 8));
}

/** Byte-swap a int16_t */
template <> int16_t byteSwap(int16_t v)
{
    uint16_t u;
    LIBERTIFF_STATIC_ASSERT(sizeof(v) == sizeof(u));
    std::memcpy(&u, &v, sizeof(u));
    u = byteSwap(u);
    std::memcpy(&v, &u, sizeof(u));
    return v;
}

/** Byte-swap a uint32_t */
template <> uint32_t byteSwap(uint32_t v)
{
    return (v >> 24) | (((v >> 16) & 0xff) << 8) | (((v >> 8) & 0xff) << 16) |
           ((v & 0xff) << 24);
}

/** Byte-swap a int32_t */
template <> int32_t byteSwap(int32_t v)
{
    uint32_t u;
    LIBERTIFF_STATIC_ASSERT(sizeof(v) == sizeof(u));
    std::memcpy(&u, &v, sizeof(u));
    u = byteSwap(u);
    std::memcpy(&v, &u, sizeof(u));
    return v;
}

/** Byte-swap a uint64_t */
template <> uint64_t byteSwap(uint64_t v)
{
    return (uint64_t(byteSwap(uint32_t(v & ~(0U)))) << 32) |
           byteSwap(uint32_t(v >> 32));
}

/** Byte-swap a int64_t */
template <> int64_t byteSwap(int64_t v)
{
    uint64_t u;
    std::memcpy(&u, &v, sizeof(u));
    u = byteSwap(u);
    std::memcpy(&v, &u, sizeof(u));
    return v;
}

/** Byte-swap a float */
template <> float byteSwap(float v)
{
    uint32_t u;
    LIBERTIFF_STATIC_ASSERT(sizeof(v) == sizeof(u));
    std::memcpy(&u, &v, sizeof(u));
    u = byteSwap(u);
    std::memcpy(&v, &u, sizeof(u));
    return v;
}

/** Byte-swap a double */
template <> double byteSwap(double v)
{
    uint64_t u;
    LIBERTIFF_STATIC_ASSERT(sizeof(v) == sizeof(u));
    std::memcpy(&u, &v, sizeof(u));
    u = byteSwap(u);
    std::memcpy(&v, &u, sizeof(u));
    return v;
}
}  // namespace LIBERTIFF_NS

namespace LIBERTIFF_NS
{
/** Interface to read from a file. */
class FileReader
{
  public:
    virtual ~FileReader() = default;

    /** Return file size in bytes */
    virtual uint64_t size() const = 0;

    /** Read 'count' bytes from offset 'offset' into 'buffer' and
     * return the number of bytes actually read.
     */
    virtual size_t read(uint64_t offset, size_t count, void *buffer) const = 0;
};
}  // namespace LIBERTIFF_NS

namespace LIBERTIFF_NS
{
/** Read context: associates a file, and the byte ordering of the TIFF file */
class ReadContext
{
  public:
    /** Constructor */
    ReadContext(const std::shared_ptr<const FileReader> &file,
                bool mustByteSwap)
        : m_file(file), m_mustByteSwap(mustByteSwap)
    {
    }

    /** Return if values of more than 1-byte must be byte swapped.
     * To be only taken into account when reading pixels. Tag values are
     * automatically byte-swapped */
    inline bool mustByteSwap() const
    {
        return m_mustByteSwap;
    }

    /** Return file size */
    inline uint64_t size() const
    {
        return m_file->size();
    }

    /** Read count raw bytes at offset into buffer */
    void read(uint64_t offset, size_t count, void *buffer, bool &ok) const
    {
        if (m_file->read(offset, count, buffer) != count)
            ok = false;
    }

    /** Read single value at offset */
    template <class T> T read(uint64_t offset, bool &ok) const
    {
#if __cplusplus >= 201703L
        static_assert(
            std::is_same_v<T, uint8_t> || std::is_same_v<T, int8_t> ||
            std::is_same_v<T, uint16_t> || std::is_same_v<T, int16_t> ||
            std::is_same_v<T, uint32_t> || std::is_same_v<T, int32_t> ||
            std::is_same_v<T, uint64_t> || std::is_same_v<T, int64_t> ||
            std::is_same_v<T, float> || std::is_same_v<T, double>);
#endif

        T res = 0;
        if (m_file->read(offset, sizeof(res), &res) != sizeof(res))
        {
            ok = false;
            return 0;
        }
        if LIBERTIFF_CONSTEXPR (sizeof(T) > 1)
        {
            if (m_mustByteSwap)
                res = byteSwap(res);
        }
        return res;
    }

    /** Read a unsigned rational (type == Type::Rational) */
    template <class T = uint32_t>
    double readRational(uint64_t offset, bool &ok) const
    {
        const auto numerator = read<T>(offset, ok);
        const auto denominator = read<T>(offset + sizeof(T), ok);
        if (denominator == 0)
        {
            ok = false;
            return std::numeric_limits<double>::quiet_NaN();
        }
        return double(numerator) / denominator;
    }

    /** Read a signed rational (type == Type::SRational) */
    double readSignedRational(uint64_t offset, bool &ok) const
    {
        return readRational<int32_t>(offset, ok);
    }

    /** Read length bytes at offset (typically for ASCII tag) as a string */
    std::string readString(std::string &res, uint64_t offset, size_t length,
                           bool &ok) const
    {
        res.resize(length);
        if (length > 0 && m_file->read(offset, length, &res[0]) != length)
        {
            ok = false;
            res.clear();
            return res;
        }
        // Strip trailing nul byte if found
        if (length > 0 && res.back() == 0)
            res.pop_back();
        return res;
    }

    /** Read length bytes at offset (typically for ASCII tag) as a string */
    std::string readString(uint64_t offset, size_t length, bool &ok) const
    {
        std::string res;
        readString(res, offset, length, ok);
        return res;
    }

    /** Read an array of count values starting at offset */
    template <class T>
    void readArray(std::vector<T> &array, uint64_t offset, size_t count,
                   bool &ok) const
    {
#if __cplusplus >= 201703L
        static_assert(
            std::is_same_v<T, uint8_t> || std::is_same_v<T, int8_t> ||
            std::is_same_v<T, uint16_t> || std::is_same_v<T, int16_t> ||
            std::is_same_v<T, uint32_t> || std::is_same_v<T, int32_t> ||
            std::is_same_v<T, uint64_t> || std::is_same_v<T, int64_t> ||
            std::is_same_v<T, float> || std::is_same_v<T, double>);
#endif

        array.resize(count);
        const size_t countBytes = count * sizeof(T);
        if (count > 0 &&
            m_file->read(offset, countBytes, &array[0]) != countBytes)
        {
            ok = false;
            array.clear();
        }
        if LIBERTIFF_CONSTEXPR (sizeof(T) > 1)
        {
            if (m_mustByteSwap)
            {
                if LIBERTIFF_CONSTEXPR (std::is_same<T, float>::value)
                {
                    uint32_t *uint32Array =
                        reinterpret_cast<uint32_t *>(array.data());
                    for (size_t i = 0; i < count; ++i)
                    {
                        uint32Array[i] = byteSwap(uint32Array[i]);
                    }
                }
                else if LIBERTIFF_CONSTEXPR (std::is_same<T, double>::value)
                {
                    uint64_t *uint64Array =
                        reinterpret_cast<uint64_t *>(array.data());
                    for (size_t i = 0; i < count; ++i)
                    {
                        uint64Array[i] = byteSwap(uint64Array[i]);
                    }
                }
                else
                {
                    for (size_t i = 0; i < count; ++i)
                    {
                        array[i] = byteSwap(array[i]);
                    }
                }
            }
        }
    }

    /** Read an array of count values starting at offset */
    template <class T>
    std::vector<T> readArray(uint64_t offset, size_t count, bool &ok) const
    {
        std::vector<T> array;
        readArray(array, offset, count, ok);
        return array;
    }

  private:
    const std::shared_ptr<const FileReader> m_file;
    const bool m_mustByteSwap;
};
}  // namespace LIBERTIFF_NS

namespace LIBERTIFF_NS
{
/** TIFF tag codes */
enum class Tag
{
    NewSubfileType = 254,
    SubfileType = 255,

    // Base line and extended TIFF tags
    ImageWidth = 256,
    ImageLength = 257,
    BitsPerSample = 258,
    Compression = 259,
    PhotometricInterpretation = 262,
    ImageDescription = 270,
    StripOffsets = 273,
    SamplesPerPixel = 277,
    RowsPerStrip = 278,
    StripByteCounts = 279,
    PlanarConfiguration = 284,
    Predictor = 317,
    TileWidth = 322,
    TileLength = 323,
    TileOffsets = 324,
    TileByteCounts = 325,
    SampleFormat = 339,

    // GeoTIFF tags
    GeoTIFFPixelScale = 33550,
    GeoTIFFTiePoints = 33922,
    GeoTIFFGeoKeyDirectory = 34735,
};

/** TIFF data types */
enum class Type
{
    Byte = 1,  /*! Unsigned 8-bit integer */
    ASCII = 2, /*! Character */
    Short = 3, /*! Unsigned 16-bit integer */
    Long = 4,  /*! Unsigned 32-bit integer */
    Rational =
        5,     /*! Positive number as a ratio of two unsigned 32-bit integers */
    SByte = 6, /*! Signed 8-bit integer */
    Undefined = 7, /*! Untyped 8-bit data */
    SShort = 8,    /*! Signed 16-bit integer */
    SLong = 9,     /*! Signed 32-bit integer */
    SRational =
        10,      /*! Signed number as a ratio of two signed 32-bit integers */
    Float = 11,  /*! 32-bit IEEE-754 floating point number */
    Double = 12, /*! 64-bit IEEE-754 floating point number */

    // BigTIFF additions
    Long8 = 16,  /*! Unsigned 64-bit integer */
    SLong8 = 17, /*! Signed 64-bit integer */
    IFD8 = 18,   /*! Unsigned 64-bit IFD offset */
};

/** Values of the PlanarConfiguration tag */
enum class PlanarConfiguration
{
    Contiguous = 1, /*! Single image plane */
    Separate = 2,   /*! Separate planes per sample */
};

/** Values of the PhotometricInterpretation tag */
enum class PhotometricInterpretation
{
    MinIsWhite = 0,
    MinIsBlack = 1,
    RGB = 2,
    Palette = 3,
    Mask = 4,
    Separated = 5,
    YCbCr = 6,
    Reserved_7 = 7,
    CIELab = 8,
    ICCLab = 9,
    ITULab = 10,
    // NOTE: if adding new values, edit processTag()
    // CFA = 32803,
    // Log2L = 32844,
    // Log2LUV = 32845,
};

/** Compression methods */

enum class Compression
{
    Unknown = 0,  //! libertiff special marker
    None = 1,
    CCITT_RLE = 2,
    CCITT_FAX3 = 3,
    CCITT_FAX4 = 4,
    LZW = 5,
    OldJPEG = 6,
    JPEG = 7,
    Deflate = 8,
    // NOTE: if adding new values, edit processTag()
};

/** Sample format */
enum class SampleFormat
{
    UnsignedInt = 1,
    SignedInt = 2,
    IEEEFP = 3,
    Void = 4,
    ComplexInt = 5,
    ComplexIEEEFP = 6
};

/** Content of a tag entry in a Image File Directory (IFD) */
struct TagEntry
{
    uint16_t tag = 0;    // of type Tag when recognized
    uint16_t type = 0;   // of type Type when recognized
    uint64_t count = 0;  // number of values in the tag

    // Inline values. Only valid if value_offset == 0.
    // The actual number in the arrays is count
    union
    {
        std::array<char, 8> charValues;
        std::array<uint8_t, 8> uint8Values;
        std::array<int8_t, 8> int8Values;
        std::array<uint16_t, 4> uint16Values;
        std::array<int16_t, 4> int16Values;
        std::array<uint32_t, 2> uint32Values;
        std::array<int32_t, 2> int32Values;
        std::array<float, 2> float32Values;
        std::array<double, 1>
            float64Values;  // Valid for Double, Rational, SRational
        std::array<uint64_t, 1> uint64Values = {0};
        std::array<int64_t, 1> int64Values;
    };

    uint64_t value_offset = 0;         // 0 for inline values
    bool invalid_value_offset = true;  // whether value_offset is invalid
};

// clang-format off
/** Return the size in bytes of a data type */
static uint32_t getSize(Type type)
{
    switch (type)
    {
        case Type::Byte:      return 1;
        case Type::ASCII:     return 1;
        case Type::Short:     return 2;
        case Type::Long:      return 4;
        case Type::Rational:  return 8;  // 2 Long
        case Type::SByte:     return 1;
        case Type::Undefined: break;
        case Type::SShort:    return 2;
        case Type::SLong:     return 4;
        case Type::SRational: return 8;  // 2 SLong
        case Type::Float:     return 4;
        case Type::Double:    return 8;
        case Type::Long8:     return 8;
        case Type::SLong8:    return 8;
        case Type::IFD8:      return 8;
    }
    return 1;
}

// clang-format on

/** Return the size in bytes of a data type, or 0 if unknown */
static uint32_t getSize(uint16_t type)
{
    if (type >= static_cast<uint16_t>(Type::Byte) &&
        type <= static_cast<uint16_t>(Type::Double))
    {
        return getSize(static_cast<Type>(type));
    }
    if (type >= static_cast<uint16_t>(Type::Long8) &&
        type <= static_cast<uint16_t>(Type::IFD8))
    {
        return 8;
    }
    return 0;
}

/** Represents a TIFF Image File Directory (IFD). */
class Image
{
  public:
    /** Constructor. Should not be called directly. Use the open() method */
    Image(const std::shared_ptr<const ReadContext> &rc, bool isBigTIFF)
        : m_rc(rc), m_isBigTIFF(isBigTIFF)
    {
    }

    /** Return read context */
    const std::shared_ptr<const ReadContext> &readContext() const
    {
        return m_rc;
    }

    /** Return whether the file is BigTIFF (if false, classic TIFF) */
    inline bool isBigTIFF() const
    {
        return m_isBigTIFF;
    }

    /** Return if values of more than 1-byte must be byte swapped.
     * To be only taken into account when reading pixels. Tag values are
     * automatically byte-swapped */
    inline bool mustByteSwap() const
    {
        return m_rc->mustByteSwap();
    }

    /** Return the offset of the next IFD (to pass to Image::open()),
         * or 0 if there is no more */
    inline uint64_t nextImageOffset() const
    {
        return m_nextImageOffset;
    }

    /** Return width of the image in pixels */
    inline uint32_t width() const
    {
        return m_width;
    }

    /** Return height of the image in pixels */
    inline uint32_t height() const
    {
        return m_height;
    }

    /** Return number of bits per sample */
    inline uint32_t bitsPerSample() const
    {
        return m_bitsPerSample;
    }

    /** Return number of samples (a.k.a. channels, bands) per pixel */
    inline uint32_t samplesPerPixel() const
    {
        return m_samplesPerPixel;
    }

    /** Return planar configuration */
    inline PlanarConfiguration planarConfiguration() const
    {
        return m_planarConfiguration;
    }

    /** Return planar configuration */
    inline PhotometricInterpretation photometricInterpretation() const
    {
        return m_photometricInterpretation;
    }

    /** Return compression method used */
    inline Compression compression() const
    {
        return m_compression;
    }

    /** Return predictor value (used for Deflate, LZW, ZStd, etc. compression) */
    inline uint32_t predictor() const
    {
        return m_predictor;
    }

    /** Return sample format */
    inline SampleFormat sampleFormat() const
    {
        return m_sampleFormat;
    }

    /** Return the number of rows per strip */
    inline uint32_t rowsPerStrip() const
    {
        return m_rowsPerStrip;
    }

    /** Return the number of strips/tiles.
     * Return 0 if inconsistent values between ByteCounts and Offsets arrays. */
    inline uint64_t strileCount() const
    {
        return m_strileCount;
    }

    /** Return whether image is tiled */
    inline bool isTiled() const
    {
        return m_isTiled;
    }

    /** Return tile width */
    inline uint32_t tileWidth() const
    {
        return m_tileWidth;
    }

    /** Return tile width */
    inline uint32_t tileHeight() const
    {
        return m_tileHeight;
    }

    /** Return number of tiles per row */
    uint32_t tilesPerRow() const
    {
        if (m_tileWidth > 0)
        {
            return uint32_t((uint64_t(m_width) + m_tileWidth - 1) /
                            m_tileWidth);
        }
        return 0;
    }

    /** Return number of tiles per column */
    uint32_t tilesPerCol() const
    {
        if (m_tileHeight > 0)
        {
            return uint32_t((uint64_t(m_height) + m_tileHeight - 1) /
                            m_tileHeight);
        }
        return 0;
    }

    /** Convert a tile coordinate (xtile, ytile, bandIdx) to a flat index */
    uint64_t tileCoordinateToIdx(uint32_t xtile, uint32_t ytile,
                                 uint32_t bandIdx, bool &ok) const
    {
        if (m_isTiled && m_tileWidth > 0 && m_tileHeight > 0)
        {
            const auto lTilesPerRow = tilesPerRow();
            const auto lTilesPerCol = tilesPerCol();
            if (xtile >= lTilesPerRow || ytile >= lTilesPerCol)
            {
                ok = false;
                return 0;
            }
            auto idx = uint64_t(ytile) * lTilesPerRow + xtile;
            if (bandIdx &&
                m_planarConfiguration == PlanarConfiguration::Separate)
            {
                idx += uint64_t(bandIdx) * lTilesPerCol * lTilesPerRow;
            }
            return idx;
        }
        ok = false;
        return 0;
    }

    /** Return the offset of strip/tile of index idx */
    uint64_t strileOffset(uint64_t idx, bool &ok) const
    {
        return readUIntTag(m_strileOffsetsTag, idx, ok);
    }

    /** Return the offset of a tile from its coordinates */
    uint64_t tileOffset(uint32_t xtile, uint32_t ytile, uint32_t bandIdx,
                        bool &ok) const
    {
        const auto idx = tileCoordinateToIdx(xtile, ytile, bandIdx, ok);
        return ok ? strileOffset(idx, ok) : 0;
    }

    /** Return the byte count of strip/tile of index idx */
    uint64_t strileByteCount(uint64_t idx, bool &ok) const
    {
        return readUIntTag(m_strileByteCountsTag, idx, ok);
    }

    /** Return the offset of a tile from its coordinates */
    uint64_t tileByteCount(uint32_t xtile, uint32_t ytile, uint32_t bandIdx,
                           bool &ok) const
    {
        const auto idx = tileCoordinateToIdx(xtile, ytile, bandIdx, ok);
        return ok ? strileByteCount(idx, ok) : 0;
    }

    /** Return the list of tags */
    inline const std::vector<TagEntry> &tags() const
    {
        return m_tags;
    }

    /** Return the (first) tag corresponding to a code, or nullptr if not found */
    const TagEntry *tag(uint16_t tagCode) const
    {
        for (const auto &tag : m_tags)
        {
            if (tag.tag == tagCode)
                return &tag;
        }
        return nullptr;
    }

    /** Return the (first) tag corresponding to a code, or nullptr if not found */
    const TagEntry *tag(Tag tagCode) const
    {
        return tag(static_cast<uint16_t>(tagCode));
    }

    /** Read an ASCII tag as a string */
    std::string readTagAsString(const TagEntry &tag, bool &ok) const
    {
        if (tag.type == static_cast<uint16_t>(Type::ASCII))
        {
            if (tag.value_offset)
            {
                return readContext()->readString(tag.value_offset, tag.count,
                                                 ok);
            }
            if (tag.count)
            {
                std::string res;
                res.resize(static_cast<size_t>(tag.count));
                std::memcpy(&res[0], tag.charValues.data(), res.size());
                if (res.back() == 0)
                    res.pop_back();
                return res;
            }
        }
        ok = false;
        return std::string();
    }

    /** Read a numeric tag as a vector */
    template <class T>
    std::vector<T> readTagAsVector(const TagEntry &tag, bool &ok) const
    {
        if (tag.value_offset)
        {
            return readContext()->readArray<T>(tag.value_offset, tag.count, ok);
        }
        if LIBERTIFF_CONSTEXPR (std::is_same<T, uint16_t>::value)
        {
            if (tag.type == static_cast<uint16_t>(Type::Short))
            {
                std::vector<T> v;
                v.insert(v.end(), tag.uint16Values.data(),
                         tag.uint16Values.data() +
                             static_cast<size_t>(tag.count));
                return v;
            }
        }
        else if LIBERTIFF_CONSTEXPR (std::is_same<T, uint32_t>::value)
        {
            if (tag.type == static_cast<uint16_t>(Type::Long))
            {
                std::vector<T> v;
                v.insert(v.end(), tag.uint32Values.data(),
                         tag.uint32Values.data() +
                             static_cast<size_t>(tag.count));
                return v;
            }
        }
        ok = false;
        return {};
    }

    /** Returns a new Image instance for the IFD starting at offset imageOffset */
    template <bool isBigTIFF>
    static std::unique_ptr<const Image>
    open(const std::shared_ptr<const ReadContext> &rc,
         const uint64_t imageOffset,
         const std::set<uint64_t> &alreadyVisitedImageOffsets =
             std::set<uint64_t>())
    {
        // To prevent infinite looping on corrupted files
        if (imageOffset == 0 || alreadyVisitedImageOffsets.find(imageOffset) !=
                                    alreadyVisitedImageOffsets.end())
        {
            return nullptr;
        }

        auto image = LIBERTIFF_NS::make_unique<Image>(rc, isBigTIFF);

        image->m_alreadyVisitedImageOffsets = alreadyVisitedImageOffsets;
        image->m_alreadyVisitedImageOffsets.insert(imageOffset);

        bool ok = true;
        int tagCount = 0;
        uint64_t offset = imageOffset;
        if LIBERTIFF_CONSTEXPR (isBigTIFF)
        {
            const auto tagCount64Bit = rc->read<uint64_t>(offset, ok);
            // Artificially limit to the same number of entries as ClassicTIFF
            if (tagCount64Bit > std::numeric_limits<uint16_t>::max())
                return nullptr;
            tagCount = static_cast<int>(tagCount64Bit);
            offset += sizeof(uint64_t);
        }
        else
        {
            tagCount = rc->read<uint16_t>(offset, ok);
            offset += sizeof(uint16_t);
        }
        if (!ok)
            return nullptr;
        image->m_tags.reserve(tagCount);
        for (int i = 0; i < tagCount; ++i)
        {
            TagEntry entry;

            // Read tag code
            entry.tag = rc->read<uint16_t>(offset, ok);
            offset += sizeof(uint16_t);

            // Read tag data type
            entry.type = rc->read<uint16_t>(offset, ok);
            offset += sizeof(uint16_t);

            // Read number of values
            if LIBERTIFF_CONSTEXPR (isBigTIFF)
            {
                auto count = rc->read<uint64_t>(offset, ok);
                entry.count = count;
                offset += sizeof(count);
            }
            else
            {
                auto count = rc->read<uint32_t>(offset, ok);
                entry.count = count;
                offset += sizeof(count);
            }

            uint32_t singleValue = 0;
            bool singleValueFitsInUInt32 = false;
            if (entry.count)
            {
                if LIBERTIFF_CONSTEXPR (isBigTIFF)
                {
                    image->ParseTagEntryDataOrOffset<uint64_t>(
                        entry, offset, singleValueFitsInUInt32, singleValue,
                        ok);
                }
                else
                {
                    image->ParseTagEntryDataOrOffset<uint32_t>(
                        entry, offset, singleValueFitsInUInt32, singleValue,
                        ok);
                }
            }
            if (!ok)
                return nullptr;

            image->processTag(entry, singleValueFitsInUInt32, singleValue);

            image->m_tags.push_back(entry);
        }

        image->finalTagProcessing();

        if LIBERTIFF_CONSTEXPR (isBigTIFF)
            image->m_nextImageOffset = rc->read<uint64_t>(offset, ok);
        else
            image->m_nextImageOffset = rc->read<uint32_t>(offset, ok);

        image->m_openFunc = open<isBigTIFF>;

        return std::unique_ptr<const Image>(image.release());
    }

    /** Returns a new Image instance at the next IFD, or nullptr if there is none */
    std::unique_ptr<const Image> next() const
    {
        return m_openFunc(m_rc, m_nextImageOffset,
                          m_alreadyVisitedImageOffsets);
    }

  private:
    const std::shared_ptr<const ReadContext> m_rc;
    std::unique_ptr<const Image> (*m_openFunc)(
        const std::shared_ptr<const ReadContext> &, const uint64_t,
        const std::set<uint64_t> &) = nullptr;

    std::set<uint64_t> m_alreadyVisitedImageOffsets{};
    uint64_t m_nextImageOffset = 0;
    uint32_t m_width = 0;
    uint32_t m_height = 0;
    uint32_t m_bitsPerSample = 0;
    uint32_t m_samplesPerPixel = 0;
    uint32_t m_rowsPerStrip = 0;
    Compression m_compression = Compression::None;
    SampleFormat m_sampleFormat = SampleFormat::UnsignedInt;
    PlanarConfiguration m_planarConfiguration = PlanarConfiguration::Contiguous;
    PhotometricInterpretation m_photometricInterpretation =
        PhotometricInterpretation::MinIsBlack;
    uint32_t m_predictor = 0;

    const bool m_isBigTIFF;
    bool m_isTiled = false;
    uint32_t m_tileWidth = 0;
    uint32_t m_tileHeight = 0;
    uint64_t m_strileCount = 0;

    std::vector<TagEntry> m_tags{};
    const TagEntry *m_strileOffsetsTag = nullptr;
    const TagEntry *m_strileByteCountsTag = nullptr;

    Image(const Image &) = delete;
    Image &operator=(const Image &) = delete;

    /** Process tag */
    void processTag(const TagEntry &entry, bool singleValueFitsInUInt32,
                    uint32_t singleValue)
    {
        if (singleValueFitsInUInt32)
        {
            switch (entry.tag)
            {
                case static_cast<uint16_t>(Tag::ImageWidth):
                    m_width = singleValue;
                    break;

                case static_cast<uint16_t>(Tag::ImageLength):
                    m_height = singleValue;
                    break;

                case static_cast<uint16_t>(Tag::Compression):
                {
                    if (singleValue >=
                            static_cast<uint32_t>(Compression::None) &&
                        singleValue <=
                            static_cast<uint32_t>(Compression::Deflate))
                    {
                        m_compression = static_cast<Compression>(singleValue);
                    }
                    else
                        m_compression = Compression::Unknown;
                    break;
                }

                case static_cast<uint16_t>(Tag::SamplesPerPixel):
                    m_samplesPerPixel = singleValue;
                    break;

                case static_cast<uint16_t>(Tag::RowsPerStrip):
                    m_rowsPerStrip = singleValue;
                    break;

                case static_cast<uint16_t>(Tag::PlanarConfiguration):
                {
                    if (singleValue >= static_cast<uint32_t>(
                                           PlanarConfiguration::Contiguous) &&
                        singleValue <= static_cast<uint32_t>(
                                           PlanarConfiguration::Separate))
                    {
                        m_planarConfiguration =
                            static_cast<PlanarConfiguration>(singleValue);
                    }
                    break;
                }

                case static_cast<uint16_t>(Tag::PhotometricInterpretation):
                {
                    if (singleValue >=
                            static_cast<uint32_t>(
                                PhotometricInterpretation::MinIsWhite) &&
                        singleValue <= static_cast<uint32_t>(
                                           PhotometricInterpretation::ITULab))
                    {
                        m_photometricInterpretation =
                            static_cast<PhotometricInterpretation>(singleValue);
                    }
                    break;
                }

                case static_cast<uint16_t>(Tag::Predictor):
                    m_predictor = singleValue;
                    break;

                case static_cast<uint16_t>(Tag::TileWidth):
                    m_tileWidth = singleValue;
                    break;

                case static_cast<uint16_t>(Tag::TileLength):
                    m_tileHeight = singleValue;
                    break;

                default:
                    break;
            }
        }

        if (entry.count && (entry.type == static_cast<uint16_t>(Type::Byte) ||
                            entry.type == static_cast<uint16_t>(Type::Short) ||
                            entry.type == static_cast<uint16_t>(Type::Long)))
        {
            // Values of those 2 tags are repeated per sample, but should be
            // at the same value.
            if (entry.tag == static_cast<uint16_t>(Tag::SampleFormat))
            {
                bool localOk = true;
                const auto sampleFormat =
                    static_cast<uint32_t>(readUIntTag(&entry, 0, localOk));
                if (localOk &&
                    sampleFormat >=
                        static_cast<uint32_t>(SampleFormat::UnsignedInt) &&
                    sampleFormat <=
                        static_cast<uint32_t>(SampleFormat::ComplexIEEEFP))
                {
                    m_sampleFormat = static_cast<SampleFormat>(sampleFormat);
                }
            }
            else if (entry.tag == static_cast<uint16_t>(Tag::BitsPerSample))
            {
                bool localOk = true;
                m_bitsPerSample =
                    static_cast<uint32_t>(readUIntTag(&entry, 0, localOk));
            }
        }
    }

    /** Final tag processing */
    void finalTagProcessing()
    {
        if ((m_strileOffsetsTag = tag(Tag::TileOffsets)) &&
            (m_strileByteCountsTag = tag(Tag::TileByteCounts)) &&
            m_strileOffsetsTag->count == m_strileByteCountsTag->count)
        {
            m_isTiled = true;
            m_strileCount = m_strileOffsetsTag->count;
        }
        else if ((m_strileOffsetsTag = tag(Tag::StripOffsets)) &&
                 (m_strileByteCountsTag = tag(Tag::StripByteCounts)) &&
                 m_strileOffsetsTag->count == m_strileByteCountsTag->count)
        {
            m_strileCount = m_strileOffsetsTag->count;
        }
    }

    /** Read a value from a byte/short/long/long8 array tag */
    uint64_t readUIntTag(const TagEntry *tag, uint64_t idx, bool &ok) const
    {
        if (tag && idx < tag->count)
        {
            if (tag->type == static_cast<uint16_t>(Type::Byte))
            {
                if (tag->count <= (m_isBigTIFF ? 8 : 4))
                {
                    return tag->uint8Values[size_t(idx)];
                }
                return m_rc->read<uint8_t>(
                    tag->value_offset + sizeof(uint8_t) * idx, ok);
            }
            else if (tag->type == static_cast<uint16_t>(Type::Short))
            {
                if (tag->count <= (m_isBigTIFF ? 4 : 2))
                {
                    return tag->uint16Values[size_t(idx)];
                }
                return m_rc->read<uint16_t>(
                    tag->value_offset + sizeof(uint16_t) * idx, ok);
            }
            else if (tag->type == static_cast<uint16_t>(Type::Long))
            {
                if (tag->count <= (m_isBigTIFF ? 2 : 1))
                {
                    return tag->uint32Values[size_t(idx)];
                }
                return m_rc->read<uint32_t>(
                    tag->value_offset + sizeof(uint32_t) * idx, ok);
            }
            else if (m_isBigTIFF &&
                     tag->type == static_cast<uint16_t>(Type::Long8))
            {
                if (tag->count <= 1)
                {
                    return tag->uint64Values[size_t(idx)];
                }
                return m_rc->read<uint64_t>(
                    tag->value_offset + sizeof(uint64_t) * idx, ok);
            }
        }
        ok = false;
        return 0;
    }

    template <class DataOrOffsetType>
    void ParseTagEntryDataOrOffset(TagEntry &entry, uint64_t &offset,
                                   bool &singleValueFitsInUInt32,
                                   uint32_t &singleValue, bool &ok)
    {
        LIBERTIFF_STATIC_ASSERT(
            (std::is_same<DataOrOffsetType, uint32_t>::value ||
             std::is_same<DataOrOffsetType, uint64_t>::value));

        const uint32_t dataTypeSize = getSize(entry.type);
        if (dataTypeSize == 0)
        {
            return;
        }

        // There are 2 cases:
        // - either the number of values for the data type can fit
        //   in the next DataOrOffsetType bytes
        // - or it cannot, and then the next DataOrOffsetType bytes are an offset
        //   to the values
        if (dataTypeSize > sizeof(DataOrOffsetType) / entry.count)
        {
            // Out-of-line values. We read a file offset
            entry.value_offset = m_rc->read<DataOrOffsetType>(offset, ok);
            const uint64_t byteCount = uint64_t(dataTypeSize) * entry.count;

            // Size of tag data beyond which we check the tag position and size
            // w.r.t the file size.
            constexpr uint32_t THRESHOLD_CHECK_FILE_SIZE = 10 * 1000 * 1000;

            entry.invalid_value_offset =
                (byteCount > THRESHOLD_CHECK_FILE_SIZE &&
                 (entry.value_offset >= m_rc->size() ||
                  entry.value_offset > m_rc->size() - byteCount));
        }
        else if (dataTypeSize == sizeof(uint8_t))
        {
            // Read up to 4 (classic) or 8 (BigTIFF) inline bytes
            m_rc->read(offset, size_t(entry.count), &entry.uint8Values[0], ok);
            if (entry.count == 1 &&
                entry.type == static_cast<uint16_t>(Type::Byte))
            {
                singleValueFitsInUInt32 = true;
                singleValue = entry.uint8Values[0];
            }
        }
        else if (dataTypeSize == sizeof(uint16_t))
        {
            // Read up to 2 (classic) or 4 (BigTIFF) inline bytes
            for (uint32_t idx = 0; idx < entry.count; ++idx)
            {
                entry.uint16Values[idx] =
                    m_rc->read<uint16_t>(offset + idx * sizeof(uint16_t), ok);
            }
            if (entry.count == 1 &&
                entry.type == static_cast<uint16_t>(Type::Short))
            {
                singleValueFitsInUInt32 = true;
                singleValue = entry.uint16Values[0];
            }
        }
        else if (dataTypeSize == sizeof(uint32_t))
        {
            // Read up to 1 (classic) or 2 (BigTIFF) inline bytes
            entry.uint32Values[0] = m_rc->read<uint32_t>(offset, ok);
            if (entry.count == 1 &&
                entry.type == static_cast<uint16_t>(Type::Long))
            {
                singleValueFitsInUInt32 = true;
                singleValue = entry.uint32Values[0];
            }
            if LIBERTIFF_CONSTEXPR (std::is_same<DataOrOffsetType,
                                                 uint64_t>::value)
            {
                if (entry.count == 2)
                {
                    entry.uint32Values[1] =
                        m_rc->read<uint32_t>(offset + sizeof(uint32_t), ok);
                }
            }
        }
        else if LIBERTIFF_CONSTEXPR (std::is_same<DataOrOffsetType,
                                                  uint64_t>::value)
        {
            if (dataTypeSize == sizeof(uint64_t))
            {
                // Read one inline 64-bit value
                if (entry.type == static_cast<uint16_t>(Type::Rational))
                    entry.float64Values[0] = m_rc->readRational(offset, ok);
                else if (entry.type == static_cast<uint16_t>(Type::SRational))
                    entry.float64Values[0] =
                        m_rc->readSignedRational(offset, ok);
                else
                    entry.uint64Values[0] = m_rc->read<uint64_t>(offset, ok);
            }
            else
            {
                assert(false);
            }
        }
        else
        {
            // fprintf(stderr, "Unexpected case: tag=%u, dataType=%u, count=%u\n", entry.tag, entry.type, entry.count);
            assert(false);
        }

        offset += sizeof(DataOrOffsetType);
    }
};

/** Open a TIFF file and return its first Image File Directory
 */
template <bool acceptBigTIFF = true>
std::unique_ptr<const Image> open(const std::shared_ptr<const FileReader> &file)
{
    unsigned char signature[2] = {0, 0};
    (void)file->read(0, 2, signature);
    const bool littleEndian = signature[0] == 'I' && signature[1] == 'I';
    const bool bigEndian = signature[0] == 'M' && signature[1] == 'M';
    if (!littleEndian && !bigEndian)
        return nullptr;

    const bool mustByteSwap = littleEndian ^ isHostLittleEndian();

    auto rc = std::make_shared<ReadContext>(file, mustByteSwap);
    bool ok = true;
    const int version = rc->read<uint16_t>(2, ok);
    constexpr int CLASSIC_TIFF_VERSION = 42;
    if (version == CLASSIC_TIFF_VERSION)
    {
        const auto firstImageOffset = rc->read<uint32_t>(4, ok);
        return Image::open<false>(rc, firstImageOffset, {});
    }
    else if LIBERTIFF_CONSTEXPR (acceptBigTIFF)
    {
        constexpr int BIGTIFF_VERSION = 43;
        if (version == BIGTIFF_VERSION)
        {
            const auto byteSizeOfOffsets = rc->read<uint16_t>(4, ok);
            if (byteSizeOfOffsets != 8)
                return nullptr;
            const auto zeroWord = rc->read<uint16_t>(6, ok);
            if (zeroWord != 0 || !ok)
                return nullptr;
            const auto firstImageOffset = rc->read<uint64_t>(8, ok);
            return Image::open<true>(rc, firstImageOffset, {});
        }
    }

    return nullptr;
}
}  // namespace LIBERTIFF_NS

#ifdef LIBERTIFF_C_FILE_READER
#include <cstdio>
#include <mutex>

namespace LIBERTIFF_NS
{
/** Interface to read from a FILE* handle */
class CFileReader final : public FileReader
{
  public:
    explicit CFileReader(FILE *f) : m_f(f)
    {
    }

    ~CFileReader() override
    {
        fclose(m_f);
    }

    uint64_t size() const override
    {
        std::lock_guard<std::mutex> oLock(m_oMutex);
        fseek(m_f, 0, SEEK_END);
        return ftell(m_f);
    }

    size_t read(uint64_t offset, size_t count, void *buffer) const override
    {
        std::lock_guard<std::mutex> oLock(m_oMutex);
        if (fseek(m_f, static_cast<long>(offset), SEEK_SET) != 0)
            return 0;
        return fread(buffer, 1, count, m_f);
    }

  private:
    FILE *const m_f;
    mutable std::mutex m_oMutex{};

    CFileReader(const CFileReader &) = delete;
    CFileReader &operator=(const CFileReader &) = delete;
};
}  // namespace LIBERTIFF_NS
#endif

#endif  // LIBERTIFF_HPP_INCLUDED
