# read in solution, which may be split into multiple files

foreach $i (6..7) {
    print "Reading in beam$i\n";
    gfx read node Cantilever$i.part0.exnode region beam$i;
    gfx read elem Cantilever$i.part0.exelem region beam$i;

    # define deformed geometry

    if ($i == 6) {
            gfx create window 1
    }

    # display deformed geometry
    gfx define faces egroup /beam$i
    gfx modify g_element "/beam$i" lines coordinate DeformedGeometry select_on material default selected_material default_selected

    gfx modify g_element "/beam$i" surfaces coordinate DeformedGeometry select_on material tissue selected_material default_selected render_shaded

    gfx modify g_element "/beam$i" node_points coordinate DeformedGeometry glyph sphere General size "2*2*2" centre 0,0,0 font default select_on material default selected_material default_selected

    # display undeformed lines
    gfx modify g_element "/beam$i" lines select_on material green selected_material default_selected

    gfx modify g_element "/" point  glyph axes general size "75*60*60" centre 0,0,0 font default select_on material default selected_material default_selected;

    gfx modify g_element "/" points domain_datapoints tessellation default_points LOCAL glyph sphere size "1*1*1" offset 0,0,0 font default select_on material default selected_material default_selected render_shaded;
    gfx modify g_element "/" points domain_datapoints tessellation default_points LOCAL glyph arrow size "1*1*1" offset 0,0,0 font default select_on material default selected_material default_selected render_shaded;

    gfx modify g_element "/" point NORMALISED_WINDOW_FIT_LEFT general size "0.5*0.5*0.2" centre -0.4,0.4,0.0 select_on material default selected_material default;
    # Removed 'glyph colour_bar' text from between '_LEFT' and 'general'

    gfx edit scene
    gfx modify window 1 set antialias 2
    gfx modify window 1 view parallel eye_point 20 -200 20 interest_point 20 20 20 up_vector 0 0 1 view_angle 40 near_clipping_plane 1.5 far_clipping_plane 700 relative_viewport ndc_placement -1 1 2 2 viewport_coordinates 0 0 1 1
}
