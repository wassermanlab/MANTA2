""" Module to handle searching the MANTA2 database for TFBS impacted be SNVs.
This module may be used in both the web application and stand alone scripts.
As such is doesn't do any DB connection management. This is presumed to be
handled outside the scope of this module.
"""

#import sys
#import bson

def search_variants(db, variants):
    """Search the MANTA2 database for TFBS impacted at any position given
    in the list variants and return a list of variant impacts.
    """

    snv_impacts = []
    errors = []
    for var in variants:
        snv_id = var['id']
        chrom = var['chrom']
        pos = var['position']
        ref_allele = var['ref_allele']
        alt_allele = var['alt_allele']

        query = {
            'chrom' : chrom,
            'snvs.pos' : pos
        }

        for tfbs_snv in db.tfbs_snvs.find(query):
            snvs = tfbs_snv['snvs']
            for snv in snvs:
                if snv['pos'] == pos:
                    if ref_allele == '.':
                        #
                        # If the input variant file format doesn't allow for
                        # the specification of the reference allele (e.g. BED
                        # format doesn't have a column to specify this), then
                        # set the reference allele to the actual reference
                        # allele stored in the MANTA DB for this SNV position.
                        #
                        ref_allele = snv['ref_allele']
                    
                    if snv['ref_allele'] != ref_allele:
                        #
                        # If the reference allele given in the input file
                        # doesn't match that in the database. Record the
                        # problem and continue with the DB reference genome
                        # assumed to be correct.
                        #
                        errors.append("WARNING: Reference allele mismatch at location chr{2}:{5} affecting transcription factor {0} ({1}) binding site chr{2}:{3}-{4}! The input SNV reference allele was given as '{6}' but the MANTA2 database records it as '{7}'. The MANTA2 database reference allele is assumed to be correct.".format(
                            tfbs_snv['jaspar_tf_name'],
                            tfbs_snv['matrix_id'],
                            chrom, tfbs_snv['start'], tfbs_snv['end'],
                            pos, ref_allele, snv['ref_allele'])
                        )

                        ref_allele = snv['ref_allele']


                    if alt_allele in snv:
                        impact = snv[alt_allele]

                        snv_impacts.append(
                            {
                                'snv_id'         : snv_id,
                                'chrom'          : chrom,
                                'position'       : pos,
                                'ref_allele'     : ref_allele,
                                'alt_allele'     : alt_allele,
                                'matrix_id'      : tfbs_snv['matrix_id'],
                                'jaspar_tf_name' : tfbs_snv['jaspar_tf_name'],
                                'start1'         : tfbs_snv['start'],
                                'end1'           : tfbs_snv['end'],
                                'strand1'        : tfbs_snv['strand'],
                                'abs_score1'     : tfbs_snv['abs_score'],
                                'rel_score1'     : tfbs_snv['rel_score'],
                                'start2'         : impact['start'],
                                'end2'           : impact['end'],
                                'strand2'        : impact['strand'],
                                'abs_score2'     : impact['abs_score'],
                                'rel_score2'     : impact['rel_score'],
                                'impact'         : impact['impact']
                            }
                        )

                        break
                    else:
                        errors.append("WARNING: Variant allele {0} at location chr{1}:{2} affecting transcription factor {3} ({4}) binding site chr{1}:{5}-{6} was NOT found - ignoring!".format(
                            alt_allele, chrom, pos,
                            tfbs_snv['jaspar_tf_name'], tfbs_snv['matrix_id'],
                            tfbs_snv['start'], tfbs_snv['end'])
                        )
                        continue

    return {'snv_impacts' : snv_impacts, 'errors' : errors}
