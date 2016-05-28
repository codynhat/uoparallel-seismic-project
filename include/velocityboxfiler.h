////////////////////////////////////////////////////////////////////////////////
// velocityboxfiler.h - 2016.05.07 - Atlee Brink
////////////////////////////////////////////////////////////////////////////////
//
// Functions for reading/writing velocity box volumes.
// Includes a function vbfileloadtext(..) to read text-format files.
//
//
// Requires:
//   point3d.h
//   velocitybox.h
//
//
// Public Datatype:
//   struct VBOXOPENFILE
//
//
// Public Functions:
//   vbfileloadtext
//
//   vbfileloadbinary 
//   vbfilestorebinary
//
//   vbfileopenbinary
//   vbfileloadbinarysubset
//   vbfileclosebinary
//
//
// Private Datatype:
//   union VBOX4BYTES
//
//
// Private Functions:
//   vbfilechecksum
//
//   vbfileread4bytes
//   vbfilewrite4bytes
//
//   vbfilereversebytes
//   vbfileread4bytesreversed
//   vbfilewrite4bytesreversed
//
//
////////////////////////////////////////////////////////////////////////////////

#pragma once


////////////////////////////////////////////////////////////////////////////////
// includes
////////////////////////////////////////////////////////////////////////////////

#include "point3d.h"
#include "velocitybox.h"


#include <stdint.h>
#include <stdio.h>
#include <string.h>


////////////////////////////////////////////////////////////////////////////////
// structs
////////////////////////////////////////////////////////////////////////////////

struct VBOXOPENFILE {
    FILE *file;    // if non-NULL, then assume file is successfully open
    long datapos;  // file position of start of float field
    int is_little_endian;
    struct POINT3D min, max, dims;
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
vbfileloadtext (
    struct VELOCITYBOX *vbox,
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
// allocates vbox->box.flat and sets everything appropriately;
// on error: returns 0
// on success: returns non-zero
{
    const char *fn = "vbfileloadtext";

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

        if( !vboxalloc( vbox, min, max ) ) {
            // couldn't allocate memory for floatbox
            fprintf (
                stderr,
                "%s: in file %s: unable to allocate memory for a VELOCITYBOX from "
                "(%d, %d, %d) to (%d, %d, %d), volume %zu\n",
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
        size_t numlines = p3dcalcvolume( vbox->box.size );
        size_t l;
        for( l = 0; l < numlines; l++ ) {
            int x, y, z;
            float vel;
            if( 4 != fscanf( infile, "%d,%d,%d,%f\n", &x, &y, &z, &vel ) ) {
                // problem parsing a line
                fprintf( stderr, "%s: I am confused by line %zu in %s\n",
                    fn, l+1, filename );
                vboxfree( vbox );
                fclose( infile );
                return 0;
            }
            vbox->box.flat[ l ] = vel;
        }
    }

    fclose( infile );

    // success
    return 1;
}


inline extern
union VBOX4BYTES
vbfilereversebytes (
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


void
vbfilechecksum (
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


int
vbfilewrite4bytes (
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

    vbfilechecksum( checksum, fb );

    // success
    return 1;
}


int
vbfilewrite4bytesreversed (
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
    union VBOX4BYTES rb = vbfilereversebytes( fb );

    if( 4 != fwrite( &rb, 1, 4, outfile ) ) {
        fprintf( stderr, "%s: error writing to file %s at position %ld\n",
            fn, filename, ftell( outfile ) );
        return 0;
    }

    vbfilechecksum( checksum, rb );

    // success
    return 1;
}


int
vbfilestorebinary (
    const char *filename,
    struct VELOCITYBOX vbox
)
// see VBOXFORMAT.txt
// on error: returns 0
// on success: returns non-zero
{
    const char *fn = "vbfilestorebinary";

    if( vbox.box.flat == NULL ) {
        fprintf( stderr, "%s: provided vbox is empty\n", fn );
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
    fb.c4[0] = (int8_t)'v';
    fb.c4[1] = (int8_t)'b';
    fb.c4[2] = (int8_t)'o';
    fb.c4[3] = (int8_t)'x';
    if( !vbfilewrite4bytes( fn, outfile, filename, fb, &checksum ) ) {
        fclose( outfile );
        return 0;
    }
    
    // total number of values to write
    size_t numvals = p3dcalcvolume( vbox.box.size );

    // detect endianness of this machine: choose path accordingly
    int err = 0;
    if( fb.u32 == 0x786f6276 ) {
        // little-endian: write values directly

        // ox, oy, oz
        fb.i32 = vbox.min.x;
        err |= !vbfilewrite4bytes( fn, outfile, filename, fb, &checksum );
        fb.i32 = vbox.min.y;
        err |= !vbfilewrite4bytes( fn, outfile, filename, fb, &checksum );
        fb.i32 = vbox.min.z;
        err |= !vbfilewrite4bytes( fn, outfile, filename, fb, &checksum );

        // check for errors
        if( err ) {
            fclose( outfile );
            return 0;
        }

        // nx, ny, nz
        fb.i32 = vbox.box.size.x;
        err |= !vbfilewrite4bytes( fn, outfile, filename, fb, &checksum );
        fb.i32 = vbox.box.size.y;
        err |= !vbfilewrite4bytes( fn, outfile, filename, fb, &checksum );
        fb.i32 = vbox.box.size.z;
        err |= !vbfilewrite4bytes( fn, outfile, filename, fb, &checksum );

        // check for errors
        if( err ) {
            fclose( outfile );
            return 0;
        }

        // flat array of velocities
        for( size_t i = 0; i < numvals; i++ ) {
            fb.f32 = vbox.box.flat[i];
            if( !vbfilewrite4bytes( fn, outfile, filename, fb, &checksum ) ) {
                fclose( outfile );
                return 0;
            }
        }

        // write final checksum value
        fb.u32 = checksum;
        if( !vbfilewrite4bytes( fn, outfile, filename, fb, &checksum ) ) {
            fclose( outfile );
            return 0;
        }

    } else {
        // big-endian: must reverse byte order before writing

        // ox, oy, oz
        fb.i32 = vbox.min.x;
        err |= !vbfilewrite4bytesreversed( fn, outfile, filename, fb, &checksum );
        fb.i32 = vbox.min.y;
        err |= !vbfilewrite4bytesreversed( fn, outfile, filename, fb, &checksum );
        fb.i32 = vbox.min.z;
        err |= !vbfilewrite4bytesreversed( fn, outfile, filename, fb, &checksum );

        // nx, ny, nz
        fb.i32 = vbox.box.size.x;
        err |= !vbfilewrite4bytesreversed( fn, outfile, filename, fb, &checksum );
        fb.i32 = vbox.box.size.y;
        err |= !vbfilewrite4bytesreversed( fn, outfile, filename, fb, &checksum );
        fb.i32 = vbox.box.size.z;
        err |= !vbfilewrite4bytesreversed( fn, outfile, filename, fb, &checksum );

        // check for errors
        if( err ) {
            fclose( outfile );
            return 0;
        }

        // flat array of velocities
        for( size_t i = 0; i < numvals; i++ ) {
            fb.f32 = vbox.box.flat[i];
            if( !vbfilewrite4bytesreversed( fn, outfile, filename, fb, &checksum ) ) {
                fclose( outfile );
                return 0;
            }
        }

        // write final checksum value
        fb.u32 = checksum;
        if( !vbfilewrite4bytesreversed( fn, outfile, filename, fb, &checksum ) ) {
            fclose( outfile );
            return 0;
        }
    }

    fclose( outfile );

    // success
    return 1;
}


int
vbfileread4bytes (
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

    vbfilechecksum( checksum, *fb );

    // success
    return 1;
}


int
vbfileread4bytesreversed (
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

    vbfilechecksum( checksum, rb );

    *fb = vbfilereversebytes( rb );

    return 1;
}


int                                 // returns: 0 on error, else non-0
vbfileopenbinary (
    struct VBOXOPENFILE *vbfile,    // out: stores metadata of open file
    const char *filename            // in: see VBOXFORMAT.txt for details
)
// note: when done with 'vbfile', be sure to run vbfileclosebinary( &vbfile )
{
    const char *fn = "vbfileopenbinary";

    vbfile->file = NULL;

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
    if( !vbfileread4bytes( fn, infile, filename, &fb, &checksum ) ) {
        fprintf( stderr, "%s: input file %s is not a vbox binary file, or is corrupted\n",
            fn, filename );
        fclose( infile );
        return 0;
    }
    if( fb.c4[0] != (int8_t)'v' || fb.c4[1] != (int8_t)'b'
    || fb.c4[2] != (int8_t)'o' || fb.c4[3] != (int8_t)'x' ) {
        // signature doesn't match
        fprintf( stderr, "%s: input file %s is not a vbox binary file, or is corrupted\n",
            fn, filename );
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
    if( (vbfile->is_little_endian = (fb.u32 == 0x786f6276)) ) { // 'vbox'
        // little-endian: read values directly

        // ox, oy, oz
        err |= !vbfileread4bytes( fn, infile, filename, &fb, &checksum );
        ox = fb.i32;
        err |= !vbfileread4bytes( fn, infile, filename, &fb, &checksum );
        oy = fb.i32;
        err |= !vbfileread4bytes( fn, infile, filename, &fb, &checksum );
        oz = fb.i32;

        // nx, ny, nz
        err |= !vbfileread4bytes( fn, infile, filename, &fb, &checksum );
        nx = fb.i32;
        err |= !vbfileread4bytes( fn, infile, filename, &fb, &checksum );
        ny = fb.i32;
        err |= !vbfileread4bytes( fn, infile, filename, &fb, &checksum );
        nz = fb.i32;
    }
    else {
        // big-endian: reverse bytes after reading

        // ox, oy, oz
        err |= !vbfileread4bytesreversed( fn, infile, filename, &fb, &checksum );
        ox = fb.i32;
        err |= !vbfileread4bytesreversed( fn, infile, filename, &fb, &checksum );
        oy = fb.i32;
        err |= !vbfileread4bytesreversed( fn, infile, filename, &fb, &checksum );
        oz = fb.i32;

        // nx, ny, nz
        err |= !vbfileread4bytesreversed( fn, infile, filename, &fb, &checksum );
        nx = fb.i32;
        err |= !vbfileread4bytesreversed( fn, infile, filename, &fb, &checksum );
        ny = fb.i32;
        err |= !vbfileread4bytesreversed( fn, infile, filename, &fb, &checksum );
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
    vbfile->file = infile;
    vbfile->datapos = ftell( infile );
    vbfile->min = p3d( ox, oy, oz );
    vbfile->dims = p3d( nx, ny, nz );
    vbfile->max = p3daddp3d( vbfile->min, p3daddval( vbfile->dims, -1 ) );
    vbfile->checksum = checksum;
    vbfile->filename = filename;

    // success
    return 1;
}


void
vbfileclosebinary (
    struct VBOXOPENFILE *vbfile
)
{
    if( vbfile == NULL || vbfile->file == NULL ) return;
    fclose( vbfile->file );
    vbfile->file = NULL;
}


int
vbfileloadbinary (
    struct VELOCITYBOX *vbox,
    const char *filename
)
// see VBOXFORMAT.txt
// on error: returns 0
// on success: returns non-zero
{
    const char *fn = "vbfileloadbinary";

    if( vbox == NULL ) return 0;

    // file metadata
    struct VBOXOPENFILE vbfile;

    // open file and read header
    if( !vbfileopenbinary( &vbfile, filename ) ) return 0;

    // prepare VELOCITYBOX to store the contents of the file
    if ( !vboxalloc( vbox, vbfile.min, vbfile.dims ) ) {
        fprintf( stderr, "%s: unable to allocate memory for a VELOCITYBOX with"
            "dimensions: %d x %d x %d\n", fn, vbfile.dims.x, vbfile.dims.y, vbfile.dims.z );
        vbfileclosebinary( &vbfile );
        return 0;
    }
    
    // this union lets us keep track of checksums and do endianness conversion
    union VBOX4BYTES fb;

    // read flat array of velocities from file and store in the VELOCITYBOX
    {
        size_t numvals = p3dcalcvolume( vbox->box.size );
        uint32_t stored_checksum;

        if( vbfile.is_little_endian ) {
            // read floats
            for( size_t i = 0; i < numvals; i++ ) {
                if( !vbfileread4bytes( fn, vbfile.file, filename, &fb, &vbfile.checksum ) ) {
                    fprintf( stderr, "%s: error reading value at byte position %zu"
                        " in %s\n", fn, ftell( vbfile.file ), filename );
                    vbfileclosebinary( &vbfile );
                    vboxfree( vbox );
                    return 0;
                }
                vbox->box.flat[i] = fb.f32;
            }
            // read stored checksum value
            uint32_t dummy_checksum;
            if( !vbfileread4bytes( fn, vbfile.file, filename, &fb, &dummy_checksum ) ) {
                fprintf( stderr, "%s: error reading stored checksum value from %s\n",
                    fn, filename );
                vbfileclosebinary( &vbfile );
                vboxfree( vbox );
                return 0;
            }
            stored_checksum = fb.u32;
        }
        else { // big endian
            // read floats
            for( size_t i = 0; i < numvals; i++ ) {
                if (
                    !vbfileread4bytesreversed (
                        fn, vbfile.file, filename, &fb, &vbfile.checksum
                    )
                ) {
                    fprintf( stderr, "%s: error reading value at byte position %zu"
                        " in %s\n", fn, ftell( vbfile.file ), filename );
                    vbfileclosebinary( &vbfile );
                    vboxfree( vbox );
                    return 0;
                }
                vbox->box.flat[i] = fb.f32;
            }
            // read stored checksum value
            uint32_t dummy_checksum;
            if( !vbfileread4bytesreversed( fn, vbfile.file, filename, &fb, &dummy_checksum ) ) {
                fprintf( stderr, "%s: error reading stored checksum value from %s\n",
                    fn, filename );
                vbfileclosebinary( &vbfile );
                vboxfree( vbox );
                return 0;
            }
            stored_checksum = fb.u32;
        }
        
        // done reading file
        vbfileclosebinary( &vbfile );

        // verify that computed vbfile.checksum matches stored checksum:
        // else something is corrupt
        if( vbfile.checksum != stored_checksum ) {
            fprintf( stderr, "%s: checksum mismatch in input file %s: suspect corruption\n",
                fn, filename );
            vboxfree( vbox );
            return 0;
        }
    }

    // success
    return 1;
}


int
vbfileloadbinarysubset (
    struct VELOCITYBOX *vbox,        // will be prepared by this function
    const struct POINT3D min,       // minimum corner of subset in file coordinates
    const struct POINT3D max,       // maximum corner of subset in file coordinates
    const struct VBOXOPENFILE vbfile // must already be open: this function doesn't close it!
)
// on error: returns 0 (either i/o error or allocation failure or something out-of-bounds)
// on success: returns non-0
// note: remember to close vbfile when you're done with it!
{
    const long valsize = sizeof(*vbox->box.flat);
    const char *fn = "vbfileloadbinarysubset";

    if( vbox == NULL ) return 0;
    if( vbfile.file == NULL ) {
        fprintf( stderr, "%s: error: source file parameter is not open\n", fn );
        return 0;
    }

    // check bounds
    if( p3disless( min, vbfile.min ) || p3dismore( max, vbfile.max ) ) {
        fprintf(
            stderr, "%s: error: file %s doesn't contain the requested subset!\n",
            fn, vbfile.filename
        );
        fprintf (
            stderr, "file: (%d, %d, %d) to (%d, %d, %d)\n",
            vbfile.min.x, vbfile.min.y, vbfile.min.z,
            vbfile.max.x, vbfile.max.y, vbfile.max.z
        );
        fprintf (
            stderr, "requested subset: (%d, %d, %d) to (%d, %d, %d)\n",
            min.x, min.y, min.z, max.x, max.y, max.z
        );
        return 0;
    }

    // prepare vbox
    if( !vboxalloc( vbox, min, max ) ) {
        fprintf(
            stderr, "%s: unable to allocate memory for a VELOCITYBOX with"
            "dimension: %d x %d x %d\n", fn, vbfile.dims.x, vbfile.dims.y, vbfile.dims.z
        );
        return 0;
    }

    // compute file strides (long because fseek uses long)
    long stridez = 1 * valsize; // if you change this, you'll need to add an fseek below
    long stridey = vbfile.dims.z * stridez;
    long stridex = vbfile.dims.y * stridey;

    // this union lets us keep track of checksums and do endianness conversion
    union VBOX4BYTES fb;

    // checksum dummy: since we're not reading the whole file, we can't use the checksum
    uint32_t dummy;

    // use different read function depending on machine endianness
    if( vbfile.is_little_endian ) {
        // little-endian: must use vbfileread4bytes(..)
        for( int x = min.x; x <= max.x; x++ ) {
            for( int y = min.y; y <= max.y; y++ ) {
                // seek to beginning of this z-strip in the file
                fseek (
                    vbfile.file,
                    vbfile.datapos +
                    (x - vbfile.min.x) * stridex +
                    (y - vbfile.min.y) * stridey +
                    (min.z - vbfile.min.z) * stridez,
                    SEEK_SET
                );
                // read one z-strip
                for( int z = min.z; z <= max.z; z++ ) {
                    if( !vbfileread4bytes( fn, vbfile.file, vbfile.filename, &fb, &dummy ) ) {
                        fprintf( stderr, "%s: error reading value at byte position %zu"
                            " in %s\n", fn, ftell( vbfile.file ), vbfile.filename );
                        vboxfree( vbox );
                        return 0;
                    }
                    vboxput( *vbox, p3d( x, y, z ), fb.f32 );
                    // assume file position was advanced exactly enough for one z
                }
            }
        }
    }
    else {
        // big-endian: must use vbfileread4bytesreversed(..)
        for( int x = min.x; x <= max.x; x++ ) {
            for( int y = min.y; y <= max.y; y++ ) {
                // seek to beginning of this z-strip in the file
                fseek (
                    vbfile.file,
                    vbfile.datapos +
                    (x - vbfile.min.x) * stridex +
                    (y - vbfile.min.y) * stridey +
                    (min.z - vbfile.min.z) * stridez,
                    SEEK_SET
                );
                // read one z-strip
                for( int z = min.z; z <= max.z; z++ ) {
                    if (
                        !vbfileread4bytesreversed (
                            fn, vbfile.file, vbfile.filename, &fb, &dummy
                        )
                    ) {
                        fprintf( stderr, "%s: error reading value at byte position %zu"
                            " in %s\n", fn, ftell( vbfile.file ), vbfile.filename );
                        vboxfree( vbox );
                        return 0;
                    }
                    vboxput( *vbox, p3d( x, y, z ), fb.f32 );
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
