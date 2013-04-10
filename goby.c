#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#ifdef HAVE_GOBY
#include <goby/C_Reads.h>
#endif
#include "goby.h"
#include "bwase.h"

#undef BWA_GOBY_WRITE_ALIGNMENT_DEBUG
#ifdef BWA_GOBY_WRITE_ALIGNMENT_DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

// The global Goby config.
goby_config_t *goby_config = NULL;

#define DEFAULT_GOBY_QUALITY_THRESHOLD_PERCENT 0.05

void goby_config_init() {
#ifdef HAVE_GOBY
    goby_config = (goby_config_t*)calloc(1, sizeof(goby_config_t));
    goby_config->output_format = OUTPUT_FORMAT_SAM;
    goby_config->output_basename = NULL;
    goby_config->input_filename = NULL;
    goby_config->input_start = 0;
    goby_config->input_end = 0;
    goby_config->which_sequence = 0;
    goby_config->num_input_files = 0;
    goby_config->num_creads_files = 0;
    goby_config->next_fafq_query_index = 0;
    goby_config->quality_threshold_percent = DEFAULT_GOBY_QUALITY_THRESHOLD_PERCENT;
#endif
}

void set_goby_config_output_basename(char *output_filename) {
#ifdef HAVE_GOBY
    if (goby_config != NULL) {
        goby_config->output_basename = strdup(output_filename);
    }
#endif
}

void bwa_print_compact_SQ(CAlignmentsWriterHelper *writerHelper, const bntseq_t *bns)
{
#ifdef HAVE_GOBY
    int i;
    for (i = 0; i < bns->n_seqs; ++i) {
        gobyAlignments_addTarget(writerHelper, i, bns->anns[i].name, bns->anns[i].len);
    }
#endif
}

#ifdef HAVE_GOBY
/**
 * Even if this match is on the reverse strand, by the time the sequences in the match
 * get here they have been reverse complemented back to the forward strand.
 * output all sequence variations (insert, delete, mutations).
 */
static void outputSequenceVariations(CAlignmentsWriterHelper *writerHelper, const char *genomic, const char *query,
        const char *quality_string, const char *fasta_query, int reverse_strand, int padding_left, int padding_right) {
  
    int nmismatches = 0, nindels = 0, i, query_i;
    char genomic_char, read_char, quality_char = '\0';
    int has_quality;
    unsigned int read_index, ref_position;
    int genomic_length = strlen(genomic);
    int fasta_length = strlen(fasta_query);
    int padded_length = padding_left + genomic_length + padding_right;
    int too_big = 0;

    // Account for alignment padding on left and right
    unsigned int *ref_positions = calloc(padded_length, sizeof(unsigned int));
    unsigned int *read_indexes = calloc(padded_length, sizeof(unsigned int));

    ref_position = 0;
    read_index = 0;
    for (i = 0; i < padded_length; i++) {
        if (i < padding_left) {
            // In alignment padding ref_position doesn't increment
            // but read position does
            read_index++;
        } else if (i >= (genomic_length + padding_left)) {
            // In alignment padding/clipping ref_position doesn't increment
            // but read position does
            read_index++;
        } else {
            genomic_char = toupper(genomic[i - padding_left]);
            if (genomic_char != '-') {
                ref_position++;
            }
            if (reverse_strand) {
                read_char = toupper(query[genomic_length - (i - padding_left) -  1]);
            } else {
                read_char = toupper(query[i - padding_left]);
            }
            if (read_char != '-') {
                read_index++;
            }
        }
        ref_positions[i] = ref_position;
        if (reverse_strand) {
            read_indexes[padded_length - i - 1] = read_index;
        } else {
            read_indexes[i] = read_index;
        }
    }

    debug(
        gobyAlignments_debugSequences(writerHelper, 1, (char *) genomic, (char *) query, padding_left, padding_right);
        fprintf(stderr, "::  pos=");
        for (i = 0; i < padded_length; i++) {
            fprintf(stderr, "%lu", ref_positions[i] % 10);
        }
        fprintf(stderr, "\n::   ri=");
        for (i = 0; i < padded_length; i++) {
            fprintf(stderr, "%lu", read_indexes[i] % 10);
        }
        fprintf(stderr, "\n");
       
        fprintf(stderr, "ref with positions\n");
        for (i = 0; i < padded_length; i++) {
            if (i < padding_left) {
                genomic_char = '_';
            } else if (i >= (genomic_length + padding_left)) {
                genomic_char = '_';
            } else {
                genomic_char = genomic[i - padding_left];
            }
            fprintf(stderr, "%02d:%c:%02lu  ", i, genomic_char, ref_positions[i]);
        }
        fprintf(stderr, "\n");
    
        fprintf(stderr, "read with positions\n");
        for (i = 0; i < padded_length; i++) {
            if (i < padding_left) {
                read_char = '_';
            } else if (i >= (padding_left + genomic_length)) {
                read_char = '_';
            } else {
                read_char = query[i - padding_left];
            }
            fprintf(stderr, "%02d:%c:%02lu  ", i, read_char, read_indexes[i]);
        }
        fprintf(stderr, "\n");
        fprintf(stderr, "padding_left=%d, padding_right=%d\n", padding_left, padding_right);
    )

    for (query_i = padding_left; query_i < padding_left + genomic_length; query_i++) {
        ref_position = ref_positions[query_i];
        read_index = read_indexes[query_i];
        i = query_i - padding_left;
        if (read_index > fasta_length) {
            too_big = 1;
        }
        genomic_char = toupper(genomic[i]);
        read_char = toupper(query[i]);
        if (quality_string != NULL && read_char != '-') {
            quality_char = quality_string[i];
            has_quality = 1;
        } else {
            has_quality = 0;
        }
        if (genomic_char != read_char) {
            gobyAlEntry_addSequenceVariation(writerHelper, read_index, ref_position, genomic_char, read_char, has_quality, quality_char);
            if (genomic_char == '-' || read_char == '-') {
                nindels++;
            } else {
                nmismatches++;
            }
        }
    }

    if (too_big) {
        fprintf(stderr, " *** read_index [%u] or ref_position [%u] is too large! ***\n",
            read_index, ref_position);
        fprintf(stderr, ">%u\n", gobyAlEntry_getQueryIndex(writerHelper));
        if (fasta_query) {
            fprintf(stderr, "%s\n", fasta_query);
        } else {
            fprintf(stderr, "fasta_query was NULL\n");
        }
    }

    gobyAlEntry_setNumberOfMismatches(writerHelper,nmismatches);
    gobyAlEntry_setScoreInt(writerHelper,genomic_length - nindels - nmismatches);

    free(ref_positions);
    free(read_indexes);

    return;
}

/**
 * Not used.
 */
static int keepEntryFilter(CSamHelper *samHelper, uint32_t queryLength) {
    if (1) {
        return 1;
    }
    long sumDifferences = samHelper->numMisMatches + samHelper->numInsertions + samHelper->numDeletions;
    double diffsAllowed = ((double) queryLength) * goby_config->quality_threshold_percent;
    // unambigious rounding function, .5 rounds up.
    diffsAllowed = floor(diffsAllowed + 0.5);
    if (sumDifferences <= diffsAllowed) {
        return 1;
    } else {
        return 0;
    }
}
#endif

/**
 * This uses the same logic as bwa_print_compact1 to determine how
 * many times the queryIndex WILL be written for the given p.
 * CURRENTLY only ambiguity==1 is supported in the logic here and in
 * bwa_print_compact1 so this returns 0 or 1.
 */
unsigned int gobyCountIndexOccurrences(const bntseq_t *bns, bwa_seq_t *p, unsigned int ambiguity) {
#ifdef HAVE_GOBY
    // If bwa_print_compact1 WILL write the entry p, this returns 1,
    // otherwise this returns 0.
    if (p->type == BWA_TYPE_NO_MATCH)  {
        // Goby doesn't store "self-unmapped"
        return 0;
    }
    int lengthRefInAlignment = pos_end(p) - p->pos;

    // get targetIndex
    int targetIndex;
    bns_coor_pac2real(bns, p->pos, lengthRefInAlignment, &targetIndex);
    if (p->pos + lengthRefInAlignment - bns->anns[targetIndex].offset > bns->anns[targetIndex].len) {
        // flag UNMAP as this alignment bridges two adjacent reference sequences
        // Goby doesn't store "self-unmapped"
        return 0;
    }
    if ((p->n_multi + 1) > ambiguity) {
        return 0;
    }
    return 1;
#else
    return 0;
#endif
}

/**
 * @param mate is only defined for paired-end matching.
 */
void bwa_print_compact1(
        CAlignmentsWriterHelper *writerHelper, const bntseq_t *bns,
        bwa_seq_t *p, const bwa_seq_t *mate, int mode,
        unsigned int queryIndexOccurrences, unsigned int ambiguity) {
#ifdef HAVE_GOBY
    // TODO do we need bwa_rg_id this for anything
    // TODO   extern char *bwa_rg_id;
    int j;
    unsigned int queryIndex;
    if (goby_config->num_creads_files) {
        queryIndex = strtoul(p->name, NULL, 10);
    } else {
        queryIndex = goby_config->next_fafq_query_index++;
    }
    gobyAlignments_observeQueryIndex(writerHelper, queryIndex);
    if (p->type == BWA_TYPE_NO_MATCH)  {
        // Goby doesn't store "self-unmapped"
        return;
    }

    int targetIndex, flag = p->extra_flag;
    int lengthRefInAlignment = pos_end(p) - p->pos;

    // get targetIndex
    bns_coor_pac2real(bns, p->pos, lengthRefInAlignment, &targetIndex);
    if (p->pos + lengthRefInAlignment - bns->anns[targetIndex].offset > bns->anns[targetIndex].len) {
        // flag UNMAP as this alignment bridges two adjacent reference sequences
        // Goby doesn't store "self-unmapped"
        return;
    }

    if (p->strand) {
        flag |= SAM_FSR; // self on the reverse strand
    }
    if (mate) {
        if (mate->type != BWA_TYPE_NO_MATCH) {
            if (mate->strand) {
                flag |= SAM_FMR;  // mate on the reverse strand
            }
        } else {
            flag |= SAM_FMU; // mate-unmapped
        }
    }

    // TODO: should p->clip_len be taken into account here or just p->full_len

    debug(fprintf(stderr,"---SAM---------------------------------\n");)
    CSamHelper *samHelper = samHelper_getResetSamHelper(writerHelper);
    samHelper_setMd(samHelper, p->md);
    samHelper_setQueryTranslate(samHelper, p->seq, p->qual, p->full_len, p->strand);
    if (p->cigar) {
        for (j = 0; j != p->n_cigar; ++j) {
            int cigarLength = __cigar_len(p->cigar[j]);
            char cigarOp = "MIDS"[__cigar_op(p->cigar[j])];
            samHelper_addCigarItem(samHelper, cigarLength, cigarOp);
        }
    } else {
        samHelper_addCigarItem(samHelper, p->full_len, 'M');
    }
    // only after reads, md, and cigar have been set, we can reconstruct the ref and reads
    // The query is now also scored, numIndels and numMisMatches set.
    samHelper_constructRefAndQuery(samHelper);

    // Determine fragment index, used with Paired alignments
    unsigned int fragmentIndex = 0;
    unsigned int m_fragmentIndex = 0;
    if (flag & SAM_FR1) {
        fragmentIndex = 0;
        m_fragmentIndex = 1;
    } else if (flag & SAM_FR2) {
        fragmentIndex = 1;
        m_fragmentIndex = 0;
    }
    // Don't +1 on the startPosition because Goby wants it 0-based
    unsigned int startPosition = (int)(p->pos - bns->anns[targetIndex].offset);

    // Mate targetIndex and startPosition, if there is a matching mate.
    int m_targetIndex = 0;
    unsigned int m_startPosition = 0;
    if (mate && mate->type != BWA_TYPE_NO_MATCH) {
        // Mate location
        bns_coor_pac2real(bns, mate->pos, mate->len, &m_targetIndex);
        m_startPosition = (int)(mate->pos - bns->anns[m_targetIndex].offset);
    }

    debug(
        fprintf(stderr, "QUERY: direction=%c\tqueryIndex=%u\tfragmentIndex=%u\ttargetIndex=%u\tposition=%u\tscore=%d\tnumInsertions=%d\tnumDeletions=%d\tnumMisMatches=%d\tcigar=%s\tmd=%s\talignedLen=%d\tref=%s\tquery=%s\torig_query=%s\tqual=%s\toorig_qual=%s\n",
            p->strand ? '-' : '+',
            queryIndex,
            fragmentIndex,
            targetIndex,
            startPosition,
            samHelper->score,
            samHelper->numInsertions,
            samHelper->numDeletions,
            samHelper->numMisMatches,
            samHelper_getCigarStr(samHelper),
            p->md,
            samHelper->alignedLength,
            samHelper_constructedRef(samHelper), 
            samHelper_constructedQuery(samHelper),
            samHelper_sourceQuery(samHelper),
            samHelper_constructedQual(samHelper),
            samHelper_sourceQual(samHelper));
        if (mate && mate->type != BWA_TYPE_NO_MATCH) {
            // Mate location
            fprintf(stderr, "MATE:  direction=%c\tqueryIndex=%u\tfragmentIndex=%u\ttargetIndex=%d\tposition=%u\n\n",
                mate->strand ? '-' : '+',
                queryIndex,
                m_fragmentIndex,
                m_targetIndex,
                m_startPosition);
        }
    )

    // TODO: Should we do something with barcode?
    // TODO: if (p->bc[0]) barcode = p->bc;

    // TODO: Better support for ambiguity lies here
    // int i;
    // if (p->n_multi) {
    //     printf("\tXA:Z:");
    //     for (i = 0; i < p->n_multi; ++i) {
    //         bwt_multi1_t *q = p->multi + i;
    //         int k;
    //         lengthRefInAlignment = pos_end_multi(q, p->len) - q->pos;
    //         bns_coor_pac2real(bns, q->pos, lengthRefInAlignment, &targetIndex);
    //         printf("%s,%c%d,", bns->anns[targetIndex].name, q->strand? '-' : '+',
    //                (int)(q->pos - bns->anns[targetIndex].offset + 1));
    //         if (q->cigar) {
    //             for (k = 0; k < q->n_cigar; ++k) {
    //                 printf("%d%c", __cigar_len(q->cigar[k]), "MIDS"[__cigar_op(q->cigar[k])]);
    //             }
    //         } else {
    //             printf("%dM", p->len);
    //         }
    //         printf(",%d;", q->gap + q->mm);
    //     }
    // }


    // TODO: This setting forces ambiguity to 1, this should be
    // TODO: changed to hande actual different ambiguity settings.
    if ((p->n_multi + 1) > ambiguity) {
        // TODO: is this right? Should provide the same output as sam-to-compact --ambiguity-threshold 1
        gobyAlEntry_appendTooManyHits(writerHelper, queryIndex, samHelper->alignedLength, p->n_multi + 1);
    }  else {
        // One match
        gobyAlignments_appendEntry(writerHelper);
        gobyAlEntry_setMultiplicity(writerHelper, 1);
        gobyAlEntry_setMatchingReverseStrand(writerHelper, p->strand);
        gobyAlEntry_setQueryIndex(writerHelper,  queryIndex);
        gobyAlEntry_setTargetIndex(writerHelper, targetIndex);
        gobyAlEntry_setQueryLength(writerHelper, p->full_len);

        const char *query = samHelper_sourceQuery(samHelper);
        const char *quals = samHelper_sourceQual(samHelper);
        gobyAlEntry_setQueryAlignedLength(writerHelper,
                samHelper->alignedLength -
                samHelper->numLeftClipped - samHelper->numRightClipped - 
                samHelper->numDeletions);
        gobyAlEntry_setTargetAlignedLength(writerHelper,
                samHelper->alignedLength -
                samHelper->numLeftClipped - samHelper->numRightClipped - 
                samHelper->numInsertions);

        if (samHelper->numLeftClipped > 0) {
            gobyAlEntry_setSoftClippedLeft(writerHelper,
                    0, samHelper->numLeftClipped, query, quals);
        }
        if (samHelper->numRightClipped > 0) {
            gobyAlEntry_setSoftClippedRight(writerHelper,
                    p->full_len - samHelper->numRightClipped,
                    samHelper->numRightClipped, query, quals);
        }
        
        gobyAlEntry_setScoreInt(writerHelper, samHelper->score);
        gobyAlEntry_setPosition(writerHelper, startPosition);
        gobyAlEntry_setMappingQuality(writerHelper, p->mapQ);
        gobyAlEntry_setNumberOfMismatches(writerHelper, samHelper->numMisMatches);
        gobyAlEntry_setNumberOfIndels(writerHelper, samHelper->numInsertions + samHelper->numDeletions);
        gobyAlEntry_setFragmentIndex(writerHelper,  fragmentIndex);
        gobyAlEntry_setPairFlags(writerHelper,  flag);

        gobyAlEntry_setAmbiguity(writerHelper, ambiguity);
        gobyAlEntry_setQueryIndexOccurrences(writerHelper, queryIndexOccurrences);
        
        if (mate != NULL) {
            if (mate->type == BWA_TYPE_NO_MATCH) {
                // Mate is unmapped
                gobyAlEntry_setPlacedUnmapped(writerHelper,
                        mate->full_len, 1 /*translateQuery*/,
                        p->strand /*reverseStrand*/,
                        mate->seq, mate->qual);
            } else {
                // Mate is mapped
                gobyAlEntry_setPairFragmentIndex(writerHelper,  m_fragmentIndex);
                gobyAlEntry_setPairTargetIndex(writerHelper, m_targetIndex);
                gobyAlEntry_setPairPosition(writerHelper, m_startPosition);
            }
        }

        // TODO: replace with new version of outputSequenceVariations()
        if ((samHelper->numInsertions + samHelper->numDeletions + samHelper->numMisMatches) > 0) {
            outputSequenceVariations(writerHelper, samHelper_constructedRef(samHelper),
                samHelper_constructedQuery(samHelper), samHelper_constructedQual(samHelper),
                samHelper_sourceQuery(samHelper), p->strand, 0, 0);
        }

        writerHelper->numberOfAlignedReads++;
    }
#endif
}

void free_goby_config() {
#ifdef HAVE_GOBY
    if (goby_config != NULL) {
        if (goby_config->output_basename) {
            free(goby_config->output_basename);
            goby_config->output_basename = NULL;
        }
        free(goby_config);
        goby_config = NULL;
    }
#endif
}
