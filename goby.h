#ifndef GOBY_H_
#define GOBY_H_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>
#include "bwtaln.h"
#include "bntseq.h"

#define OUTPUT_FORMAT_SAM  0
#define OUTPUT_FORMAT_GOBY 1

#ifdef HAVE_GOBY
#include <goby/C_Reads.h> 	/* has gobyReads_..., goby_shutdownProtobuf,  */
#include <goby/C_Alignments.h> 	/* has gobyReads_..., goby_shutdownProtobuf,  */
#include <goby/C_CompactHelpers.h>
#endif

#undef BWA_GOBY_READ_SEQ_DEBUG
#ifdef BWA_GOBY_READ_SEQ_DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

typedef struct {
    int output_format;
    char *output_basename;
	char *input_filename;
	unsigned long input_start;
	unsigned long input_end;
	int which_sequence;
	double quality_threshold_percent;
    int num_input_files;
    int num_creads_files;
    unsigned int next_fafq_query_index;
} goby_config_t;

extern goby_config_t *goby_config;

#ifdef HAVE_GOBY
#else
typedef struct {
} CReadsHelper;
typedef struct {
} CAlignmentsWriterHelper;
typedef struct {
} CSamHelper;
#endif

void goby_config_init();
void set_goby_config_output_basename(char *output_filename);
void bwa_print_compact_SQ(CAlignmentsWriterHelper *writerHelper, const bntseq_t *bns);
unsigned int gobyCountIndexOccurrences(const bntseq_t *bns, bwa_seq_t *p, unsigned int ambiguity);
void bwa_print_compact1(CAlignmentsWriterHelper *writerHelper, const bntseq_t *bns,
        bwa_seq_t *p, const bwa_seq_t *mate, int mode,
        unsigned int queryIndexOccurrences, unsigned int ambiguity);
void free_goby_config();

#ifdef HAVE_GOBY

#define __KSEQ_READ_GOBY                                                \
    static int clip_warn = 0;                                           \
    static int kseq_read_compact(kseq_t *seq)                           \
    {                                                                   \
        if (!gobyReads_hasNext(seq->f->creads_helper)) {                \
            return -1;                                                  \
        }                                                               \
        char *identifier;                                               \
        char *description;                                              \
        char *sequence;                                                 \
        int sequenceLength;                                             \
        char *quality;                                                  \
        int qualityLength;                                              \
        char *psequence;                                                \
        int psequenceLength;                                            \
        char *pquality;                                                 \
        int pqualityLength;                                             \
        unsigned int readIndex = gobyReads_nextSequencePair(           \
            seq->f->creads_helper,                                      \
            &identifier, &description,                                  \
            &sequence, &sequenceLength, &quality, &qualityLength,       \
            &psequence, &psequenceLength, &pquality, &pqualityLength);  \
        if (seq->f->which_creads_sequence != 0) {                       \
            if (psequenceLength == 0) {                                 \
                fprintf(stderr, "%s %s\n",                              \
                    "ERROR! User specified to read pair sequence,",     \
                    "but length is 0. Using non-pair sequence.");       \
            } else {                                                    \
                /* Use the pair instead of the primary */               \
                sequence = psequence;                                   \
                sequenceLength = psequenceLength;                       \
                quality = pquality;                                     \
                qualityLength = pqualityLength;                         \
            }                                                           \
        }                                                               \
        if (qualityLength != 0 && (sequenceLength != qualityLength)) {  \
            if (clip_warn == 0) {                                       \
                fprintf(stderr,                                         \
                    "This warning will only appear once:\n");           \
                fprintf(stderr,                                         \
                    "!Quality and sequence length don't match %d, %d\n",\
                    sequenceLength, qualityLength);                     \
            }                                                           \
            if (qualityLength > sequenceLength) {                       \
                if (clip_warn == 0) {                                   \
                    fprintf(stderr,                                     \
                        "!Clipping qual to sequence length\n");         \
                    clip_warn = 1;                                      \
                }                                                       \
                qualityLength = sequenceLength;                         \
                quality[qualityLength] = '\0';                          \
            } else {                                                    \
                return -2;                                              \
            }                                                           \
        }                                                               \
 	seq->comment.l = seq->seq.l = seq->qual.l = seq->last_char = 0;	\
        /* Name */                                                      \
        int readIndexLength = ((readIndex == 0) ? 1 :                   \
            ((((int)log10(readIndex)) + 1)));                           \
        if (readIndexLength + 1 > seq->name.m) {                        \
            seq->name.m = readIndexLength + 1;                          \
            seq->name.s = (char*)realloc(seq->name.s, seq->name.m);     \
        }                                                               \
        sprintf(seq->name.s, "%u", readIndex);                         \
        seq->name.l = readIndexLength;                                  \
        /* Comment */                                                   \
		int descriptionLength = strlen(description);                    \
        if (descriptionLength + 1 > seq->comment.m) {                   \
            seq->comment.m = descriptionLength + 1;                     \
            seq->comment.s = (char*)realloc(                            \
                seq->comment.s, seq->comment.m);                        \
        }                                                               \
        strcpy(seq->comment.s, description);                            \
        seq->comment.l = descriptionLength;                             \
        /* Sequence */                                                  \
        if (sequenceLength + 1 > seq->seq.m) {                          \
            seq->seq.m = sequenceLength + 1;                            \
            seq->seq.s = (char*)realloc(seq->seq.s, seq->seq.m);        \
        }                                                               \
        strcpy(seq->seq.s, sequence);                                   \
        seq->seq.l = sequenceLength;                                    \
        /* Quality */                                                   \
        if (qualityLength > 0) {                                        \
            if (qualityLength + 1 > seq->qual.m) {                      \
                seq->qual.m = qualityLength + 1;                        \
                seq->qual.s = (char*)realloc(seq->qual.s, seq->qual.m); \
            }                                                           \
            strcpy(seq->qual.s, quality);                               \
            seq->qual.l = qualityLength;                                \
        }                                                               \
        return seq->seq.l;                                              \
	}

#define __OPEN_GOBY_COMPACT                                             \
        int firstChar = ks_getc(s->f);                                  \
        /* Unget that first char from the buffer. */                    \
        s->f->begin--;                                                  \
        if (firstChar == -1 && goby_config != NULL &&                   \
                goby_config->input_filename != NULL) {                  \
            /** Close original file, open as Compact Reads. */          \
            gzclose(s->f->f);                                           \
            s->f->f = NULL;                                             \
            s->f->is_creads = 1;                                        \
            goby_config->num_creads_files += 1;                         \
            s->f->which_creads_sequence = goby_config->which_sequence;  \
            fprintf(stderr,                                             \
                "Opening Goby compact-reads %s [%lu-%lu] Sequence %d\n",\
                goby_config->input_filename,                            \
                goby_config->input_start, goby_config->input_end,       \
                s->f->which_creads_sequence);                           \
            gobyReads_openReadsReaderSingleWindowed(                    \
                goby_config->input_filename,                            \
                goby_config->input_start, goby_config->input_end,       \
                &s->f->creads_helper);                                  \
            /* Convert from phred to illumina scores when reading. */   \
            gobyReads_setQualityAdjustment(s->f->creads_helper, 64);    \
        }

#else /* HAVE_GOBY */

#define __KSEQ_READ_GOBY                                                \
    static int kseq_read_compact(kseq_t *seq)                           \
    {                                                                   \
       return -1;                                                       \
    }

#define __OPEN_GOBY_COMPACT

#endif /* HAVE_GOBY */
#endif /* GOBY_H_ */
