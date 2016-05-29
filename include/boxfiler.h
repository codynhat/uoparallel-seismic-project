////////////////////////////////////////////////////////////////////////////////
// boxfiler.h - 2016.05.07 - Atlee Brink
////////////////////////////////////////////////////////////////////////////////
//
// Functions for reading/writing float box volumes.
// Includes a function boxfileloadtext(..) to read text-format files.
//
//
// Requires:
//   floatbox.h
//   point3d.h
//
//
// Public Datatype:
//   struct BOXOPENFILE
//
//
// Public Functions:
//   boxfileloadtext
//
//   boxfileloadbinary 
//   boxfilestorebinary
//
//   boxfileopenbinary
//   boxfileloadbinarysubset
//   boxfileclosebinary
//
//
// Private Datatype:
//   union VBOX4BYTES
//
//
// Private Functions:
//   boxfilechecksum
//
//   boxfileread4bytes
//   boxfilewrite4bytes
//
//   boxfilereversebytes
//   boxfileread4bytesreversed
//   boxfilewrite4bytesreversed
//
//
////////////////////////////////////////////////////////////////////////////////

#pragma once


////////////////////////////////////////////////////////////////////////////////
// includes
////////////////////////////////////////////////////////////////////////////////

#include "floatbox.h"
#include "point3d.h"


#include <stdint.h>
#include <stdio.h>
#include <string.h>


////////////////////////////////////////////////////////////////////////////////
// structs
////////////////////////////////////////////////////////////////////////////////

struct BOXOPENFILE {
    FILE *file;    // if non-NULL, then assume file is successfully open
    long datapos;  // file position of start of float field
    int is_little_endian;
    struct POINT3D min, max, size;
    uint32_t checksum;
    const char *filename;
};


////////////////////////////////////////////////////////////////////////////////
// unions
////////////////////////////////////////////////////////////////////////////////

union VBOX4BYTES {
    int8_t c4[4];
    int32_t i32;
    uint32_t u32;
    float f32; 
};


////////////////////////////////////////////////////////////////////////////////
// functions
////////////////////////////////////////////////////////////////////////////////

int
boxfileloadtext (
    struct FLOATBOX *box,
    const char *filename
)
// reads text-format velocity file of the format:
//   x_1,y_1,z_1,float_1
//   x_1,y_1,z_2,float_2
//   ...
//   x_nx,y_ny,z_nz,float_nxnynz
// example:
//   1,1,1,0.29762
//   ...
//   241,241,51,0.1771
// this function prepares 'box' and fills it with the entire volume
// on error: returns 0
// on success: returns non-zero
{
    const char *fn = "boxfileloadtext";

    FILE *infile = fopen( filename, "r" );
    if( infile == NULL ) {
        // there was a problem opening the given file
        fprintf( stderr, "%s: error opening file %s\n", fn, filename );
        return 0;
    }

    // will receive coordinates of entire volume in file
    struct POINT3D min, max;

    // read first line to get (x,y,z) of first velocity
    if( 3 != fscanf( infile, "%d,%d,%d", &min.x, &min.y, &min.z ) ) {
        // error
        fprintf( stderr, "%s: error reading first line from file %s\n",
            fn, filename );
        fclose( infile );
        return 0;
    }

    // read last line to get (x,y,z) of last velocity;
    // this is used with the first coordinates to determine box size
    {
        int minlinelen = strlen( "1,1,1,0.0" );

        if( 0 != fseek( infile, -minlinelen, SEEK_END ) ) {
            // couldn't seek to estimated last line of file: maybe file is too short
            fprintf( stderr, "%s: error seeking to estimated last line in file %s\n",
                fn, filename );
            fclose( infile );
            return 0;
        }

        // scan backward for a newline character
        for(;;) {
            int c = fgetc( infile );
            if( c == '\n' || c == '\r' ) {
                // found a newline: the file cursor is now at the start of the last line
                break;
            }
            if( c == EOF || 0 != fseek( infile, -2, SEEK_CUR ) ) {
                // some i/o error occurred
                fprintf( stderr, "%s: error scanning for last line in file %s\n",
                    fn, filename );
                fclose( infile );
                return 0;
            }
        }

        // read coordinates (assume they are correct)
        if( 3 != fscanf( infile, "%d,%d,%d", &max.x, &max.y, &max.z ) ) {
            // error
            fprintf( stderr, "%s: error reading last line from file %s\n",
                fn, filename );
            fclose( infile );
            return 0;
        }

        // check that the volume coordinates make sense
        if( p3dismore( min, max ) ) {
            fprintf( stderr, "%s: nonsense coordinates in file %s\n",
                fn, filename );
            fclose( infile );
            return 0;
        }

        if( !boxalloc( box, min, max, min, max ) ) {
            // couldn't allocate memory for floatbox
            fprintf (
                stderr,
                "%s: in file %s: unable to allocate memory for a FLOATBOX from "
                "(%d, %d, %d) to (%d, %d, %d), volume %ld\n",
                fn, filename,
                min.x, min.y, min.z, max.x, max.y, max.z,
                p3dcalcvolume( p3dsizeofregion( min, max ) )
            );
            fclose( infile );
            return 0;
        }
    }

    // starting from the beginning of the file, read and remember all velocity values
    {
        fseek( infile, 0, SEEK_SET );
        long numlines = p3dcalcvolume( box->size );
        long l;
        for( l = 0; l < numlines; l++ ) {
            int x, y, z;
            float vel;
            if( 4 != fscanf( infile, "%d,%d,%d,%f\n", &x, &y, &z, &vel ) ) {
                // problem parsing a line
                fprintf( stderr, "%s: I am confused by line %ld in %s\n",
                    fn, l+1, filename );
                boxfree( box );
                fclose( infile );
                return 0;
            }
            box->flat[ l ] = vel;
        }
    }

    fclose( infile );

    // success
    return 1;
}


union VBOX4BYTES
boxfilesigtofb (
    const char sig[4] // such as {'v','b','o','x'}
)
{
    union VBOX4BYTES fb;
    fb.c4[0] = (int8_t)sig[0];
    fb.c4[1] = (int8_t)sig[1];
    fb.c4[2] = (int8_t)sig[2];
    fb.c4[3] = (int8_t)sig[3];
    return fb;
}


uint32_t
boxfilesigval (
    const char sig[4] // such as {'v','b','o','x'}
)
{
    return
        ((uint32_t)sig[3] & 0xff) << 24 |
        ((uint32_t)sig[2] & 0xff) << 16 |
        ((uint32_t)sig[1] & 0xff) << 8 |
        ((uint32_t)sig[0] & 0xff)
    ;
}


inline extern
union VBOX4BYTES
boxfilereversebytes (
    const union VBOX4BYTES in
)
// reverses the byte-order of four bytes
{
    union VBOX4BYTES out;

    out.c4[0] = in.c4[3];
    out.c4[1] = in.c4[2];
    out.c4[2] = in.c4[1];
    out.c4[3] = in.c4[0];

    return out;
}


inline extern
void
boxfilechecksum (
    uint32_t *checksum,
    union VBOX4BYTES fb
)
// updates a checksum based on the bytes in fb;
// assumes little-endian byte order, so use this AFTER conversion, if necessary;
{
    *checksum += (uint32_t)fb.c4[0]
        + ((uint32_t)fb.c4[1] << 8)
        + ((uint32_t)fb.c4[2] << 16)
        + ((uint32_t)fb.c4[3] << 24);
}


inline extern
int
boxfilewrite4bytes (
    const char *fn,
    FILE *outfile,
    const char *filename,
    union VBOX4BYTES fb,
    uint32_t *checksum
)
// writes exactly 4 bytes to the given output FILE*;
// updates byte position counter and updates checksum
// on error: returns 0
// on success: returns non-zero
{
    if( 4 != fwrite( &fb, 1, 4, outfile ) ) {
        fprintf( stderr, "%s: error writing to file %s at position %ld\n",
            fn, filename, ftell( outfile ) );
        return 0;
    }

    boxfilechecksum( checksum, fb );

    // success
    return 1;
}


inline extern
int
boxfilewrite4bytesreversed (
    const char *fn,
    FILE *outfile,
    const char *filename,
    union VBOX4BYTES fb,
    uint32_t *checksum
)
// writes exactly 4 bytes to the given output FILE*, but in reverse byte order;
// updates byte position counter and updates checksum
// on error: returns 0
// on success: returns non-zero
{
    union VBOX4BYTES rb = boxfilereversebytes( fb );

    if( 4 != fwrite( &rb, 1, 4, outfile ) ) {
        fprintf( stderr, "%s: error writing to file %s at position %ld\n",
            fn, filename, ftell( outfile ) );
        return 0;
    }

    boxfilechecksum( checksum, rb );

    // success
    return 1;
}


int
boxfilestorebinary (
    const char *filename,
    const char signature[4], // example: {'v','b','o','x'} for a velocity box
    const struct FLOATBOX box
)
// see VBOXFORMAT.txt
// on error: returns 0
// on success: returns non-zero
{
    const char *fn = "boxfilestorebinary";

    if( box.flat == NULL ) {
        fprintf( stderr, "%s: provided box is empty\n", fn );
        return 0;
    }

    FILE *outfile = fopen( filename, "wb" );
    if( outfile == NULL ) {
        // there was a problem opening/creating the given file
        fprintf( stderr, "%s: error creating file %s\n", fn, filename );
        return 0;
    }

    // this union lets us do stuff like checksums and endianness conversion
    union VBOX4BYTES fb;

    // keep track of write position and checksum
    uint32_t checksum = 0;

    // write signature
    fb = boxfilesigtofb( signature );
    if( !boxfilewrite4bytes( fn, outfile, filename, fb, &checksum ) ) {
        fclose( outfile );
        return 0;
    }
    
    // total number of values to write
    long numvals = p3dcalcvolume( box.size );

    // detect endianness of this machine: choose path accordingly
    uint32_t sigval = boxfilesigval( signature );
    int err = 0;
    if( fb.u32 == sigval ) {
        // little-endian: write values directly

        // ox, oy, oz
        fb.i32 = box.omin.x;
        err |= !boxfilewrite4bytes( fn, outfile, filename, fb, &checksum );
        fb.i32 = box.omin.y;
        err |= !boxfilewrite4bytes( fn, outfile, filename, fb, &checksum );
        fb.i32 = box.omin.z;
        err |= !boxfilewrite4bytes( fn, outfile, filename, fb, &checksum );

        // check for errors
        if( err ) {
            fclose( outfile );
            return 0;
        }

        // nx, ny, nz
        fb.i32 = box.size.x;
        err |= !boxfilewrite4bytes( fn, outfile, filename, fb, &checksum );
        fb.i32 = box.size.y;
        err |= !boxfilewrite4bytes( fn, outfile, filename, fb, &checksum );
        fb.i32 = box.size.z;
        err |= !boxfilewrite4bytes( fn, outfile, filename, fb, &checksum );

        // check for errors
        if( err ) {
            fclose( outfile );
            return 0;
        }

        // flat array of velocities
        for( long i = 0; i < numvals; i++ ) {
            fb.f32 = box.flat[i];
            if( !boxfilewrite4bytes( fn, outfile, filename, fb, &checksum ) ) {
                fclose( outfile );
                return 0;
            }
        }

        // write final checksum value
        fb.u32 = checksum;
        if( !boxfilewrite4bytes( fn, outfile, filename, fb, &checksum ) ) {
            fclose( outfile );
            return 0;
        }
    }
    else {
        // big-endian: must reverse byte order before writing

        // ox, oy, oz
        fb.i32 = box.omin.x;
        err |= !boxfilewrite4bytesreversed( fn, outfile, filename, fb, &checksum );
        fb.i32 = box.omin.y;
        err |= !boxfilewrite4bytesreversed( fn, outfile, filename, fb, &checksum );
        fb.i32 = box.omin.z;
        err |= !boxfilewrite4bytesreversed( fn, outfile, filename, fb, &checksum );

        // nx, ny, nz
        fb.i32 = box.size.x;
        err |= !boxfilewrite4bytesreversed( fn, outfile, filename, fb, &checksum );
        fb.i32 = box.size.y;
        err |= !boxfilewrite4bytesreversed( fn, outfile, filename, fb, &checksum );
        fb.i32 = box.size.z;
        err |= !boxfilewrite4bytesreversed( fn, outfile, filename, fb, &checksum );

        // check for errors
        if( err ) {
            fclose( outfile );
            return 0;
        }

        // flat array of velocities
        for( long i = 0; i < numvals; i++ ) {
            fb.f32 = box.flat[i];
            if( !boxfilewrite4bytesreversed( fn, outfile, filename, fb, &checksum ) ) {
                fclose( outfile );
                return 0;
            }
        }

        // write final checksum value
        fb.u32 = checksum;
        if( !boxfilewrite4bytesreversed( fn, outfile, filename, fb, &checksum ) ) {
            fclose( outfile );
            return 0;
        }
    }

    fclose( outfile );

    // success
    return 1;
}


inline extern
int
boxfileread4bytes (
    const char *fn,
    FILE *infile,
    const char *filename,
    union VBOX4BYTES *fb,
    uint32_t *checksum
)
// reads exactly 4 bytes from the given input FILE*;
// updates byte position counter and updates checksum
// on error: returns 0
// on success: returns non-zero
{
    if( 4 != fread( fb, 1, 4, infile ) ) {
        fprintf( stderr, "%s: error reading from file %s at position %ld\n",
            fn, filename, ftell( infile ) );
        return 0;
    }

    boxfilechecksum( checksum, *fb );

    // success
    return 1;
}


inline extern
int
boxfileread4bytesreversed (
    const char *fn,
    FILE *infile,
    const char *filename,
    union VBOX4BYTES *fb,
    uint32_t *checksum
)
// reads exactly 4 bytes from the given input FILE*, but in reverse byte order;
// updates byte position counter and updates checksum
// on error: returns 0
// on success: returns non-zero
{
    union VBOX4BYTES rb;

    if( 4 != fread( &rb, 1, 4, infile ) ) {
        fprintf( stderr, "%s: error reading from file %s at position %ld\n",
            fn, filename, ftell( infile ) );
        return 0;
    }

    boxfilechecksum( checksum, rb );

    *fb = boxfilereversebytes( rb );

    return 1;
}


int                                 // returns: 0 on error, else non-0
boxfileopenbinary (
    struct BOXOPENFILE *boxfile,    // out: stores metadata of open file
    const char *filename,           // in: see VBOXFORMAT.txt for details
    const char signature[4]         // example: {'v','b','o','x'} for a velocity box
)
// note: when done with 'boxfile', be sure to run boxfileclosebinary( &boxfile )
{
    const char *fn = "boxfileopenbinary";

    boxfile->file = NULL;

    FILE *infile = fopen( filename, "rb" );
    if( infile == NULL ) {
        fprintf( stderr, "%s: error opening file %s\n", fn, filename );
        return 0;
    }

    // this union lets us keep track of checksums and do endianness conversion
    union VBOX4BYTES fb;

    // keep track of read position and checksum
    uint32_t checksum = 0;

    // read signature
    if( !boxfileread4bytes( fn, infile, filename, &fb, &checksum ) ) {
        fprintf( stderr, "%s: input file %s is not a box file, or is corrupted\n",
            fn, filename );
        fclose( infile );
        return 0;
    }
    union VBOX4BYTES fbsigexp = boxfilesigtofb( signature );
    if( fb.u32 != fbsigexp.u32 ) {
        // file signature doesn't match expected signature
        fprintf (
            stderr,
            "%s: input file %s doesn't have the correct signature, and may be corrupt:\n"
            "\t    file signature: \"%c%c%c%c\"\n"
            "\texpected signature: \"%c%c%c%c\"\n",
            fn, filename,
            fb.c4[0], fb.c4[1], fb.c4[2], fb.c4[3],
            fbsigexp.c4[0], fbsigexp.c4[1], fbsigexp.c4[2], fbsigexp.c4[3]
        );
        fclose( infile );
        return 0;
    }

    // global minimum coordinates of velocity volume
    int ox, oy, oz;

    // dimensions of velocity volume
    int nx, ny, nz;

    // check after reading header: if err != 0, then there was an error
    int err = 0;

    // detect endianness of this machine: choose path accordingly
    if( (boxfile->is_little_endian = (fb.u32 == boxfilesigval( signature ))) ) { // 'vbox'
        // little-endian: read values directly

        // ox, oy, oz
        err |= !boxfileread4bytes( fn, infile, filename, &fb, &checksum );
        ox = fb.i32;
        err |= !boxfileread4bytes( fn, infile, filename, &fb, &checksum );
        oy = fb.i32;
        err |= !boxfileread4bytes( fn, infile, filename, &fb, &checksum );
        oz = fb.i32;

        // nx, ny, nz
        err |= !boxfileread4bytes( fn, infile, filename, &fb, &checksum );
        nx = fb.i32;
        err |= !boxfileread4bytes( fn, infile, filename, &fb, &checksum );
        ny = fb.i32;
        err |= !boxfileread4bytes( fn, infile, filename, &fb, &checksum );
        nz = fb.i32;
    }
    else {
        // big-endian: reverse bytes after reading

        // ox, oy, oz
        err |= !boxfileread4bytesreversed( fn, infile, filename, &fb, &checksum );
        ox = fb.i32;
        err |= !boxfileread4bytesreversed( fn, infile, filename, &fb, &checksum );
        oy = fb.i32;
        err |= !boxfileread4bytesreversed( fn, infile, filename, &fb, &checksum );
        oz = fb.i32;

        // nx, ny, nz
        err |= !boxfileread4bytesreversed( fn, infile, filename, &fb, &checksum );
        nx = fb.i32;
        err |= !boxfileread4bytesreversed( fn, infile, filename, &fb, &checksum );
        ny = fb.i32;
        err |= !boxfileread4bytesreversed( fn, infile, filename, &fb, &checksum );
        nz = fb.i32;
    }

    // check for errors
    if( err ) {
        fprintf( stderr, "%s: error reading header in %s: suspect corruption\n",
            fn, filename );
        fclose( infile );
        return 0;
    }

    // as appropriate, store metadata from header
    boxfile->file = infile;
    boxfile->datapos = ftell( infile );
    boxfile->min = p3d( ox, oy, oz );
    boxfile->size = p3d( nx, ny, nz );
    boxfile->max = p3daddp3d( boxfile->min, p3daddval( boxfile->size, -1 ) );
    boxfile->checksum = checksum;
    boxfile->filename = filename;

    // success
    return 1;
}


void
boxfileclosebinary (
    struct BOXOPENFILE *boxfile
)
{
    if( boxfile == NULL || boxfile->file == NULL ) return;
    fclose( boxfile->file );
    boxfile->file = NULL;
}


int
boxfileloadbinary (
    struct FLOATBOX *box,
    const char signature[4], // example: {'v','b','o','x'} for a velocity box
    const char *filename
)
// see VBOXFORMAT.txt
// on error: returns 0
// on success: returns non-zero
{
    const char *fn = "boxfileloadbinary";

    if( box == NULL ) {
        fprintf( stderr, "%s: error: received a null-pointer instead of a FLOATBOX*\n",
            fn );
        return 0;
    }

    // file metadata
    struct BOXOPENFILE boxfile;

    // open file and read header
    if( !boxfileopenbinary( &boxfile, signature, filename ) ) return 0;

    // prepare FLOATBOX to store the contents of the file
    if ( !boxalloc( box, boxfile.min, boxfile.max, boxfile.min, boxfile.max ) ) {
        fprintf( stderr, "%s: unable to allocate memory for a FLOATBOX with"
            "dimensions: %d x %d x %d\n", fn, boxfile.size.x, boxfile.size.y, boxfile.size.z );
        boxfileclosebinary( &boxfile );
        return 0;
    }
    
    // this union lets us keep track of checksums and do endianness conversion
    union VBOX4BYTES fb;

    // read flat array of velocities from file and store in the FLOATBOX
    {
        long numvals = p3dcalcvolume( box->size );
        uint32_t stored_checksum;

        if( boxfile.is_little_endian ) {
            // read floats
            for( long i = 0; i < numvals; i++ ) {
                if( !boxfileread4bytes( fn, boxfile.file, filename, &fb, &boxfile.checksum ) ) {
                    fprintf( stderr, "%s: error reading value at byte position %ld"
                        " in %s\n", fn, ftell( boxfile.file ), filename );
                    boxfileclosebinary( &boxfile );
                    boxfree( box );
                    return 0;
                }
                box->flat[i] = fb.f32;
            }
            // read stored checksum value
            uint32_t dummy_checksum;
            if( !boxfileread4bytes( fn, boxfile.file, filename, &fb, &dummy_checksum ) ) {
                fprintf( stderr, "%s: error reading stored checksum value from %s\n",
                    fn, filename );
                boxfileclosebinary( &boxfile );
                boxfree( box );
                return 0;
            }
            stored_checksum = fb.u32;
        }
        else { // big endian
            // read floats
            for( long i = 0; i < numvals; i++ ) {
                if (
                    !boxfileread4bytesreversed (
                        fn, boxfile.file, filename, &fb, &boxfile.checksum
                    )
                ) {
                    fprintf( stderr, "%s: error reading value at byte position %ld"
                        " in %s\n", fn, ftell( boxfile.file ), filename );
                    boxfileclosebinary( &boxfile );
                    boxfree( box );
                    return 0;
                }
                box->flat[i] = fb.f32;
            }
            // read stored checksum value
            uint32_t dummy_checksum;
            if( !boxfileread4bytesreversed( fn, boxfile.file, filename, &fb, &dummy_checksum ) ) {
                fprintf( stderr, "%s: error reading stored checksum value from %s\n",
                    fn, filename );
                boxfileclosebinary( &boxfile );
                boxfree( box );
                return 0;
            }
            stored_checksum = fb.u32;
        }
        
        // done reading file
        boxfileclosebinary( &boxfile );

        // verify that computed boxfile.checksum matches stored checksum:
        // else something is corrupt
        if( boxfile.checksum != stored_checksum ) {
            fprintf( stderr, "%s: checksum mismatch in input file %s: suspect corruption\n",
                fn, filename );
            boxfree( box );
            return 0;
        }
    }

    // success
    return 1;
}


int
boxfileloadbinarysubset (
    struct FLOATBOX *box,            // will be prepared by this function
    const struct POINT3D omin,       // outer minimum in file coordinates
    const struct POINT3D omax,       // outer maximum in file coordinates
    const struct POINT3D imin,       // inner minimum in file coordinates
    const struct POINT3D imax,       // inner maximum in file coordinates
    const struct BOXOPENFILE boxfile // must already be open: this function doesn't close it!
)
// on error: returns 0 (either i/o error or allocation failure or something out-of-bounds)
// on success: returns non-0
// note: remember to close boxfile when you're done with it!
{
    const long valsize = sizeof(*box->flat);
    const char *fn = "boxfileloadbinarysubset";

    if( box == NULL ) {
        fprintf( stderr, "%s: error: received a null-pointer instead of a FLOATBOX*\n",
            fn );
        return 0;
    }
    if( boxfile.file == NULL ) {
        fprintf( stderr, "%s: error: source file parameter is not open\n", fn );
        return 0;
    }

    // check bounds
    if( p3disless( omin, boxfile.min ) || p3dismore( omax, boxfile.max ) ) {
        fprintf(
            stderr, "%s: error: file %s doesn't contain the requested subset!\n",
            fn, boxfile.filename
        );
        fprintf (
            stderr, "file: (%d, %d, %d) to (%d, %d, %d)\n",
            boxfile.min.x, boxfile.min.y, boxfile.min.z,
            boxfile.max.x, boxfile.max.y, boxfile.max.z
        );
        fprintf (
            stderr, "requested subset: (%d, %d, %d) to (%d, %d, %d)\n",
            omin.x, omin.y, omin.z, omax.x, omax.y, omax.z
        );
        return 0;
    }
    if( p3disless( imin, omin ) || p3dismore( imax, omax ) ) {
        fprintf (
            stderr, "%s: error: either imin < omin or imax > omax: should be contained!\n",
            fn
        );
        return 0;
    }

    // prepare box
    if( !boxalloc( box, omin, omax, imin, imax ) ) {
        fprintf(
            stderr, "%s: unable to allocate memory for a FLOATBOX with"
            "dimension: %d x %d x %d\n", fn, boxfile.size.x, boxfile.size.y, boxfile.size.z
        );
        return 0;
    }

    // compute file strides (long because fseek uses long)
    long stridez = 1 * valsize; // if you change this, you'll need to add an fseek below
    long stridey = boxfile.size.z * stridez;
    long stridex = boxfile.size.y * stridey;

    // this union lets us keep track of checksums and do endianness conversion
    union VBOX4BYTES fb;

    // checksum dummy: since we're not reading the whole file, we can't use the checksum
    uint32_t dummy;

    // use different read function depending on machine endianness
    if( boxfile.is_little_endian ) {
        // little-endian: must use boxfileread4bytes(..)
        for( int x = omin.x; x <= omax.x; x++ ) {
            for( int y = omin.y; y <= omax.y; y++ ) {
                // seek to beginning of this z-strip in the file
                fseek (
                    boxfile.file,
                    boxfile.datapos +
                    (x - boxfile.min.x) * stridex +
                    (y - boxfile.min.y) * stridey +
                    (omin.z - boxfile.min.z) * stridez,
                    SEEK_SET
                );
                // read one z-strip
                for( int z = omin.z; z <= omax.z; z++ ) {
                    if( !boxfileread4bytes( fn, boxfile.file, boxfile.filename, &fb, &dummy ) ) {
                        fprintf( stderr, "%s: error reading value at byte position %ld"
                            " in %s\n", fn, ftell( boxfile.file ), boxfile.filename );
                        boxfree( box );
                        return 0;
                    }
                    boxputglobal( *box, p3d( x, y, z ), fb.f32 );
                    // assume file position was advanced exactly enough for one z
                }
            }
        }
    }
    else {
        // big-endian: must use boxfileread4bytesreversed(..)
        for( int x = omin.x; x <= omax.x; x++ ) {
            for( int y = omin.y; y <= omax.y; y++ ) {
                // seek to beginning of this z-strip in the file
                fseek (
                    boxfile.file,
                    boxfile.datapos +
                    (x - boxfile.min.x) * stridex +
                    (y - boxfile.min.y) * stridey +
                    (omin.z - boxfile.min.z) * stridez,
                    SEEK_SET
                );
                // read one z-strip
                for( int z = omin.z; z <= omax.z; z++ ) {
                    if (
                        !boxfileread4bytesreversed (
                            fn, boxfile.file, boxfile.filename, &fb, &dummy
                        )
                    ) {
                        fprintf( stderr, "%s: error reading value at byte position %ld"
                            " in %s\n", fn, ftell( boxfile.file ), boxfile.filename );
                        boxfree( box );
                        return 0;
                    }
                    boxputglobal( *box, p3d( x, y, z ), fb.f32 );
                    // assume file position was advanced exactly enough for one z
                }
            }
        }
    }

    // success
    return 1;
}


////////////////////////////////////////////////////////////////////////////////
// END
////////////////////////////////////////////////////////////////////////////////
