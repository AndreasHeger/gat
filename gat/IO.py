import sys, re, os, glob

import gat
import gat.IOTools as IOTools
import gat.Experiment as E

def buildSegments( options ):
    '''load segments, annotations and workspace from parameters
    defined in *options*.

    The workspace will be split by isochores.

    returns segments, annotations and workspace.
    '''

    def expandGlobs( infiles ):
        return IOTools.flatten( [ glob.glob( x ) for x in infiles ] )
        
    options.segment_files = expandGlobs( options.segment_files )
    options.annotation_files = expandGlobs( options.annotation_files )
    options.workspace_files = expandGlobs( options.workspace_files )
    options.sample_files = expandGlobs( options.sample_files )

    ##################################################
    # arguments sanity check
    if not options.segment_files:
        raise ValueError("please specify at least one segment file" )
    if not options.annotation_files:
        raise ValueError("please specify at least one annotation file" )
    if not options.workspace_files:
        raise ValueError("please specify at least one workspace file" )

    ##################################################
    ##################################################
    ##################################################
    # process input
    def dumpStats( coll, section ):
        if section in options.output_stats or \
                "all" in options.output_stats or \
                len( [ x for x in options.output_stats if re.search( x, section ) ] ) > 0:
            coll.outputStats( E.openOutputFile( section ) )

    def dumpBed( coll, section ):
        if section in options.output_bed or \
                "all" in options.output_bed or \
                len( [ x for x in options.output_bed if re.search( x, section ) ] ) > 0:
            coll.save( E.openOutputFile( section + ".bed" ) )

    def readSegmentList( label, filenames, enable_split_tracks = False ):
        # read one or more segment files
        results = gat.IntervalCollection( name = label )
        E.info( "%s: reading tracks from %i files" % (label, len(filenames)))
        results.load( filenames, split_tracks = enable_split_tracks )
        E.info( "%s: read %i tracks from %i files" % (label, len(results), len(filenames)))
        dumpStats( results, "stats_%s_raw" % label )
        results.normalize()
        dumpStats( results, "stats_%s_normed" % label )
        return results

    # read one or more segment files
    segments = readSegmentList( "segments", options.segment_files)
    if options.ignore_segment_tracks:
        segments.merge( delete = True)
        E.info( "merged all segments into one track with %i segments" % len(segments))

    if len(segments) > 1000: 
        raise ValueError( "too many (%i) segment files - use track definitions or --ignore-segment-tracks" % len(segments) )
    
    annotations = readSegmentList( "annotations", options.annotation_files, options.enable_split_tracks )
    workspaces = readSegmentList( "workspaces", options.workspace_files, options.enable_split_tracks )

    # intersect workspaces to build a single workspace
    E.info( "collapsing workspaces" )
    dumpStats( workspaces, "stats_workspaces_input" )
    workspaces.collapse()
    dumpStats( workspaces, "stats_workspaces_collapsed" )

    # use merged workspace only, discard others
    workspaces.restrict("collapsed")

    # build isochores or intersect annotations/segments with workspace
    if options.isochore_files:
        
        # read one or more isochore files
        isochores = gat.IntervalCollection( name = "isochores" )
        E.info( "%s: reading isochores from %i files" % ("isochores", len(options.isochore_files)))
        isochores.load( options.isochore_files )
        dumpStats( isochores, "stats_isochores_raw" )

        # merge isochores and check if consistent (fully normalized)
        isochores.sort()

        # check that there are no overlapping segments within isochores
        isochores.check()

        # TODO: flag is_normalized not properly set
        isochores.normalize()

        # check that there are no overlapping segments between isochores

        # truncate isochores to workspace
        # crucial if isochores are larger than workspace.
        isochores.intersect( workspaces["collapsed"] )

    else:
        isochores = None
        
    return segments, annotations, workspaces, isochores

def applyIsochores( segments, annotations, workspaces, isochores = None ):
    '''apply isochores to segments, annotations.

    Segments and annotations are filtered to keep only those overlapping
    the workspace.

    If no isochores are given, isochores are not applied.

    returns workspace divided into isochores.
    '''

    if isochores:
        # intersect isochores and workspaces, segments and annotations
        E.info( "adding isochores to workspace" )
        workspaces.toIsochores( isochores )
        annotations.toIsochores( isochores )
        segments.toIsochores( isochores )
        
        if workspaces.sum() == 0:
            raise ValueError( "isochores and workspaces do not overlap" )
        if annotations.sum() == 0:
            raise ValueError( "isochores and annotations do not overlap" )
        if segments.sum() == 0:
            raise ValueError( "isochores and segments do not overlap" )

        dumpStats( workspaces, "stats_workspaces_isochores" )
        dumpStats( annotations, "stats_annotations_isochores" )
        dumpStats( segments, "stats_segments_isochores" )
    
        dumpBed( workspaces, "workspaces_isochores" )
        dumpBed( annotations, "annotations_isochores" )
        dumpBed( segments, "segments_isochores" )

    else:
        # intersect workspace and segments/annotations
        annotations.filter( workspaces["collapsed"] )
        segments.filter( workspaces["collapsed"] )
        
        dumpStats( annotations, "stats_annotations_pruned" )
        dumpStats( segments, "stats_segments_pruned" )

    workspace = workspaces["collapsed"] 

    if options.restrict_workspace:
        E.info( "restricting workspace" )
        # this is very cumbersome - refactor merge and collapse
        # to return an IntervalDictionary instead of adding it
        # to the list of tracks
        for x in (segments, annotations):
            if "merged" in segments:
                workspace.filter( segments["merged"] )
            else:
                segments.merge()
                workspace.filter( segments["merged"] )
                del segments[merged]

        dumpStats( workspaces, "stats_workspaces_restricted" )
        
    # segments.dump( open("segments_dump.bed", "w" ) )
    # workspaces.dump( open("workspaces_dump.bed", "w" ) )

    # output overlap stats
    # output segment densities per workspace
    for track in segments.tracks:
        workspaces.outputOverlapStats( E.openOutputFile( "overlap_%s" % track), 
                                       segments[track] )

    return workspace
