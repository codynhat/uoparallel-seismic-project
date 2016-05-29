////////////////////////////////////////////////////////////////////////////////
// example_velocityboxfiler.c - Atlee Brink
// Demonstrates and tests velocityboxfiler.h.
////////////////////////////////////////////////////////////////////////////////

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

    struct FLOATBOX box_orig;

    // load a text-format velocity file
    printf( "loading old velocity model %s...", text_velocity_file );
    fflush( stdout );

    if( !vbfileloadtext( &box_orig, text_velocity_file ) ) {
        // handle error
    }

    printf( " done.\n" );
    fflush( stdout );

    // print the old vbox metadata
    boxfprint( stdout, "old: ", "\t", box_orig );

    // don't free box_orig yet: we're going to use it below


////////////////////////////////////////////////////////////////////////////////
// writing vbox-format file
////////////////////////////////////////////////////////////////////////////////

    puts("");

    // Write box_orig to a new file using the vbox file format
    printf( "storing velocity model to vbox-format file %s...", vbox_output_file );
    fflush( stdout );

    if( !vbfilestorebinary( vbox_output_file, box_orig ) ) {
        // handle error
    }
    
    printf( " done.\n" );
    fflush( stdout );

    // don't free box_orig yet: we're going to use it below


////////////////////////////////////////////////////////////////////////////////
// reading whole vbox-format file
////////////////////////////////////////////////////////////////////////////////

    puts("");

    struct FLOATBOX box_new;

    // load the vbox file we just created
    printf( "loading velocity model from vbox-format file %s...", vbox_output_file);
    fflush( stdout );
    
    if( !vbfileloadbinary( &box_new, vbox_output_file ) ) {
        // handle error
    }

    printf( " done.\n" );
    fflush( stdout );

    // print the new vbox metadata
    boxfprint( stdout, "new: ", "\t", box_new );

    // don't free box_new yet: we're going to use it below


////////////////////////////////////////////////////////////////////////////////
// using global coordinates with boxgetglobal to compare box_orig and box_new
////////////////////////////////////////////////////////////////////////////////

    puts("");

    {
        printf( "comparing values in %s to values in %s...\n",
            text_velocity_file, vbox_output_file );
        fflush( stdout );
        int bad = 0;
        size_t count = 0;
        for( int x = box_orig.omin.x; !bad && x <= box_orig.omax.x; x++ ) {
            for( int y = box_orig.omin.y; !bad && y <= box_orig.omax.y; y++ ) {
                for( int z = box_orig.omin.z; z <= box_orig.omax.z; z++ ) {
                    count++;
                    struct POINT3D pt = {x, y, z};
                    float val_orig = boxgetglobal( box_orig, pt );
                    float val_new = boxgetglobal( box_new, pt );
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
                p3dcalcvolume( p3dsizeofregion( box_orig.omin, box_orig.omax ) ) );
        }
    }

    // free box_new, but don't free box_orig yet: we're going to use it below
    boxfree( &box_new );


////////////////////////////////////////////////////////////////////////////////
// reading subset of vbox-format file
////////////////////////////////////////////////////////////////////////////////

    puts("");

    struct VBOXOPENFILE vbfile;
    struct FLOATBOX box_subset;

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
    struct POINT3D third = { vbfile.size.x / 3, vbfile.size.y / 3, vbfile.size.z / 3 };
    struct POINT3D min = p3daddp3d( vbfile.min, third );
    struct POINT3D max = p3daddp3d( min, third );

    printf( "subset min: (%d, %d, %d)\n", min.x, min.y, min.z );
    printf( "subset max: (%d, %d, %d)\n", max.x, max.y, max.z );

    printf( "loading subset..." );
    fflush( stdout );

    // load a subset of the file into a VELOCITYBOX
    if( !vbfileloadbinarysubset( &box_subset, min, max, min, max, vbfile ) ) {
        // handle error
    }

    printf( " done.\n" );
    fflush( stdout );

    // close the open vbfile
    vbfileclosebinary( &vbfile );

    // don't free box_subset yet: we're going to use it below


////////////////////////////////////////////////////////////////////////////////
// using global coordinates to compare box_subset with box_orig
////////////////////////////////////////////////////////////////////////////////

    puts("");

    {
        printf( "comparing subset of new file to original text file...\n" );
        fflush( stdout );
        int bad = 0;
        size_t count = 0;
        for( int x = box_subset.omin.x; !bad && x <= box_subset.omax.x; x++ ) {
            for( int y = box_subset.omin.y; !bad && y <= box_subset.omax.y; y++ ) {
                for( int z = box_subset.omin.z; z <= box_subset.omax.z; z++ ) {
                    count++;
                    struct POINT3D pt = {x, y, z};
                    float val_orig = boxgetglobal( box_orig, pt );
                    float val_subset = boxgetglobal( box_subset, pt );
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
                p3dcalcvolume( p3dsizeofregion( box_subset.omin, box_subset.omax ) ) );
        }
    }

    // free heap memory used by the original and subset
    boxfree( &box_orig );
    boxfree( &box_subset );


////////////////////////////////////////////////////////////////////////////////
// done
////////////////////////////////////////////////////////////////////////////////

    puts("");

    printf( "all done.\n" );

    return 0;
}


////////////////////////////////////////////////////////////////////////////////
// END
////////////////////////////////////////////////////////////////////////////////
