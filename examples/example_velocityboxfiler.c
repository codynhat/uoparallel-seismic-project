// example_velocityboxfiler.c - Atlee Brink
// Demonstrates and tests velocityboxfiler.h.

#include "velocityboxfiler.h"

#include <stdio.h> /* for stdout, not necessary otherwise */

const char *text_velocity_file = "../docs/velocity-241-241-51.txt";
const char *vbox_output_file = "velocities.vbox";
const struct POINT3D pos = {101, 105, 25};

int main()
{
////////////////////////////////////////////////////////////////////////////////
// reading text-format file
////////////////////////////////////////////////////////////////////////////////

    struct VELOCITYBOX vbox_orig;

    // load a text-format velocity file
    printf( "loading old velocity model %s...", text_velocity_file );
    fflush( stdout );

    if( !vbfileloadtext( &vbox_orig, text_velocity_file ) ) {
        // handle error
    }

    printf( " done.\n" );
    fflush( stdout );

    // print the old vbox metadata
    vboxfprint( stdout, "old: ", "\t", vbox_orig );

    // don't free vbox_orig yet: we're going to use it below


////////////////////////////////////////////////////////////////////////////////
// writing vbox-format file
////////////////////////////////////////////////////////////////////////////////

    puts("");

    // Write vbox_orig to a new file using the vbox file format
    printf( "storing velocity model to vbox-format file %s...", vbox_output_file );
    fflush( stdout );

    if( !vbfilestorebinary( vbox_output_file, vbox_orig ) ) {
        // handle error
    }
    
    printf( " done.\n" );
    fflush( stdout );

    // don't free vbox_orig yet: we're going to use it below


////////////////////////////////////////////////////////////////////////////////
// reading whole vbox-format file
////////////////////////////////////////////////////////////////////////////////

    puts("");

    struct VELOCITYBOX vbox_new;

    // load the vbox file we just created
    printf( "loading velocity model from vbox-format file %s...", vbox_output_file);
    fflush( stdout );
    
    if( !vbfileloadbinary( &vbox_new, vbox_output_file ) ) {
        // handle error
    }

    printf( " done.\n" );
    fflush( stdout );

    // print the new vbox metadata
    vboxfprint( stdout, "new: ", "\t", vbox_new );

    // don't free vbox_new yet: we're going to use it below


////////////////////////////////////////////////////////////////////////////////
// using global coordinates with vboxget to compare vbox_orig and vbox_new
////////////////////////////////////////////////////////////////////////////////

    puts("");

    {
        printf( "comparing values in %s to values in %s...\n",
            text_velocity_file, vbox_output_file );
        fflush( stdout );
        int bad = 0;
        size_t count = 0;
        for( int x = vbox_orig.min.x; !bad && x <= vbox_orig.max.x; x++ ) {
            for( int y = vbox_orig.min.y; !bad && y <= vbox_orig.max.y; y++ ) {
                for( int z = vbox_orig.min.z; z <= vbox_orig.max.z; z++ ) {
                    count++;
                    struct POINT3D pt = {x, y, z};
                    float val_orig = vboxget( vbox_orig, pt );
                    float val_new = vboxget( vbox_new, pt );
                    if( val_orig != val_new ) {
                        printf( "value mismatch at global coordinate (%d, %d, %d)\n",
                            x, y, z );
                        printf( "%s says: %g\n", text_velocity_file, val_orig );
                        printf( "%s says: %g\n", vbox_output_file, val_new );
                        printf( "aborting!\n" );
                        bad = 1;
                        break;
                    }
                }
            }
        }
        if( !bad ) {
            printf( "values agree for %zu of %zu elements\n", count,
                p3dcalcvolume( p3dsizeofregion( vbox_orig.min, vbox_orig.max ) ) );
        }
    }

    // free vbox_new, but don't free vbox_orig yet: we're going to use it below
    vboxfree( &vbox_new );


////////////////////////////////////////////////////////////////////////////////
// reading subset of vbox-format file
////////////////////////////////////////////////////////////////////////////////

    puts("");

    struct VBOXOPENFILE vbfile;
    struct VELOCITYBOX vbox_subset;

    // open (don't read whole file yet)
    printf( "opening just-written vbox-format velocity model %s...", vbox_output_file);
    fflush( stdout );
    
    if( !vbfileopenbinary( &vbfile, vbox_output_file ) ) {
        printf( "error opening %s!\n", vbox_output_file );
        return 0;
    }
    
    printf( " done.\n" );
    fflush( stdout );

    // subset region
    struct POINT3D third = { vbfile.dims.x / 3, vbfile.dims.y / 3, vbfile.dims.z / 3 };
    struct POINT3D min = p3daddp3d( vbfile.min, third );
    struct POINT3D max = p3daddp3d( min, third );

    printf( "subset min: (%d, %d, %d)\n", min.x, min.y, min.z );
    printf( "subset max: (%d, %d, %d)\n", max.x, max.y, max.z );

    printf( "loading subset..." );
    fflush( stdout );

    // load a subset of the file into a VELOCITYBOX
    if( !vbfileloadbinarysubset( &vbox_subset, min, max, vbfile ) ) {
        // handle error
    }

    printf( " done.\n" );
    fflush( stdout );

    // close the open vbfile
    vbfileclosebinary( &vbfile );

    // don't free vbox_subset yet: we're going to use it below


////////////////////////////////////////////////////////////////////////////////
// using global coordinates to compare vbox_subset with vbox_orig
////////////////////////////////////////////////////////////////////////////////

    puts("");

    {
        printf( "comparing subset of new file to original text file...\n" );
        fflush( stdout );
        int bad = 0;
        size_t count = 0;
        for( int x = vbox_subset.min.x; !bad && x <= vbox_subset.max.x; x++ ) {
            for( int y = vbox_subset.min.y; !bad && y <= vbox_subset.max.y; y++ ) {
                for( int z = vbox_subset.min.z; z <= vbox_subset.max.z; z++ ) {
                    count++;
                    struct POINT3D pt = {x, y, z};
                    float val_orig = vboxget( vbox_orig, pt );
                    float val_subset = vboxget( vbox_subset, pt );
                    if( val_orig != val_subset ) {
                        printf( "value mismatch at global coordinate (%d, %d, %d)\n",
                            x, y, z );
                        printf( "original says: %g\n", val_orig );
                        printf( "subset says: %g\n", val_subset );
                        printf( "aborting!\n" );
                        bad = 1;
                        break;
                    }
                }
            }
        }
        if( !bad ) {
            printf( "values agree for %zu of %zu elements\n", count,
                p3dcalcvolume( p3dsizeofregion( vbox_subset.min, vbox_subset.max ) ) );
        }
    }

    // free heap memory used by the original and subset
    vboxfree( &vbox_orig );
    vboxfree( &vbox_subset );


////////////////////////////////////////////////////////////////////////////////
// done
////////////////////////////////////////////////////////////////////////////////

    puts("");

    printf( "all done.\n" );

    return 0;
}
