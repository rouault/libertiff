// Define LIBERTIFF_C_FILE_READER to be able to use libertiff::CFileReader
#define LIBERTIFF_C_FILE_READER
#include "libertiff.hpp"

#include <cinttypes>

static void dumpIFD(const libertiff::Image *tiff)
{
    printf("Width: %u\n", tiff->width());
    printf("Height: %u\n", tiff->height());
    printf("SamplesPerPixel: %u\n", tiff->samplesPerPixel());
    printf("PlanarConfiguration: %u\n",
           static_cast<uint32_t>(tiff->planarConfiguration()));
    printf("Compression: %u\n", static_cast<uint32_t>(tiff->compression()));
    printf("PhotometricInterpretation: %u\n",
           static_cast<uint32_t>(tiff->photometricInterpretation()));
    printf("Predictor: %u\n", tiff->predictor());
    printf("SampleFormat: %u\n", static_cast<uint32_t>(tiff->sampleFormat()));
    if (!tiff->isTiled())
    {
        printf("Organization: strip\n");
        printf("Number of strips: %" PRIu64 "\n", tiff->strileCount());
        printf("RowsPerStrip: %u\n", tiff->rowsPerStrip());
        printf("StripOffsets: ");
        for (uint64_t i = 0; i < tiff->strileCount(); ++i)
        {
            bool ok = true;
            if (i > 0)
                printf(",");
            printf("%" PRIu64, tiff->strileOffset(i, ok));
        }
        printf("\n");
        printf("StripByteCounts: ");
        for (uint64_t i = 0; i < tiff->strileCount(); ++i)
        {
            bool ok = true;
            if (i > 0)
                printf(",");
            printf("%" PRIu64, tiff->strileByteCount(i, ok));
        }
        printf("\n");
    }
    else
    {
        printf("Organization: tile\n");
        printf("Number of tiles: %" PRIu64 "\n", tiff->strileCount());
        printf("TileWidth: %u\n", tiff->tileWidth());
        printf("TileHeight: %u\n", tiff->tileHeight());
        printf("TilesPerRow: %u\n", tiff->tilesPerRow());
        printf("TilesPerCol: %u\n", tiff->tilesPerCol());
        printf("TileOffsets: ");
        for (uint64_t i = 0; i < tiff->strileCount(); ++i)
        {
            bool ok = true;
            if (i > 0)
                printf(",");
            printf("%" PRIu64, tiff->strileOffset(i, ok));
        }
        printf("\n");
        printf("TileByteCounts: ");
        for (uint64_t i = 0; i < tiff->strileCount(); ++i)
        {
            bool ok = true;
            if (i > 0)
                printf(",");
            printf("%" PRIu64, tiff->strileByteCount(i, ok));
        }
        printf("\n");
    }
    const auto &tags = tiff->tags();
    printf("All tags (%u):\n", static_cast<uint32_t>(tags.size()));
    for (const auto &tag : tags)
    {
        const char *typeName = "";
        if (tag.type == static_cast<uint16_t>(libertiff::Type::ASCII))
            typeName = "ASCII";
        else if (tag.type == static_cast<uint16_t>(libertiff::Type::Short))
            typeName = "SHORT";
        else if (tag.type == static_cast<uint16_t>(libertiff::Type::Long))
            typeName = "LONG";
        else if (tag.type == static_cast<uint16_t>(libertiff::Type::Long8))
            typeName = "LONG8";
        printf("  Code=%d(0x%X), Type=%s(%d), Count=%" PRIu64, tag.tag, tag.tag,
               typeName, tag.type, tag.count);
        if (!tag.invalid_value_offset)
        {
            printf(", Offset of value=%" PRIu64, tag.value_offset);
        }
        if (tag.count < 1000)
        {
            if (tag.type == static_cast<uint16_t>(libertiff::Type::ASCII))
            {
                bool ok = true;
                printf(", Value=\"%s\"",
                       tiff->readTagAsString(tag, ok).c_str());
            }
            else if (tag.type == static_cast<uint16_t>(libertiff::Type::Short))
            {
                printf(", Value%s=", tag.count > 1 ? "s" : "");
                bool ok = true;
                const auto values = tiff->readTagAsVector<uint16_t>(tag, ok);
                for (size_t i = 0; i < values.size(); ++i)
                {
                    if (i > 0)
                        printf(",");
                    printf("%u", values[i]);
                }
            }
            else if (tag.type == static_cast<uint16_t>(libertiff::Type::Long))
            {
                printf(", Value%s=", tag.count > 1 ? "s" : "");
                bool ok = true;
                const auto values = tiff->readTagAsVector<uint32_t>(tag, ok);
                for (size_t i = 0; i < values.size(); ++i)
                {
                    if (i > 0)
                        printf(",");
                    printf("%u", values[i]);
                }
            }
            else if (tag.type == static_cast<uint16_t>(libertiff::Type::Long8))
            {
                printf(", Value%s=", tag.count > 1 ? "s" : "");
                bool ok = true;
                const auto values = tiff->readTagAsVector<uint64_t>(tag, ok);
                for (size_t i = 0; i < values.size(); ++i)
                {
                    if (i > 0)
                        printf(",");
                    printf("%" PRIu64, values[i]);
                }
            }
        }

        printf("\n");
    }
}

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        fprintf(stderr, "Usage: demo my.tiff\n");
        exit(1);
    }

    FILE *f = fopen(argv[1], "rb");
    if (!f)
    {
        fprintf(stderr, "Cannot open %s\n", argv[1]);
        exit(1);
    }

    auto tiff = libertiff::open(std::make_shared<libertiff::CFileReader>(f));
    if (!tiff)
    {
        fprintf(stderr, "Not a valid TIFF file\n");
        exit(1);
    }

    printf("File: %s\n", argv[1]);
    uint32_t ifd = 0;
    do
    {
        printf("IFD %u:\n", ifd);
        dumpIFD(tiff.get());
        tiff = tiff->next();
    } while (tiff);

    return 0;
}
