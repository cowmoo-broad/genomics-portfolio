task STARAlignGenome {
    File reference_fa
    File reference_gtf

    command {
        mkdir out
        /opt/STAR-2.7.10a/bin/Linux_x86_64/STAR --runMode genomeGenerate \
            --genomeDir out \
            --genomeFastaFiles ${reference_fa} \
            --sjdbGTFfile ${reference_gtf} \
            --sjdbOverhang 62 \
    }
    output {
        Array[File] out = glob("out/*.*")
    }

    runtime {
        docker: "cowmoo/genomics-container:latest"
    }
}

task STARAlignTranscripts {
    Array[File] reference_outs
    File fastq

    command <<<
        mkdir out
        for x in ~{sep=' ' reference_outs}
           ln-s "$\{x}" out/"$\{x}"
        done;

        /opt/STAR-2.7.10a/bin/Linux_x86_64/STAR --runThreadsN 4 --genomeDir out --readFilesIn ${fastq}
    >>>

    runtime {
        docker: "cowmoo/genomics-container:latest"
    }
}

task fastqc {

    File fastq

    command {
       fastqc ${fastq}
    }

    output {
        Array[File] out = glob("*.zip")
    }

    runtime {
        docker: "cowmoo/genomics-container:latest"
    }
}

task RSEMPrepareReference {

    File exon_gtf
    File reference_fa

    command {
        export PATH=/bowtie-1.3.1-linux-x86_64/:$PATH
        rsem-prepare-reference --bowtie --gtf ${exon_gtf} ${reference_fa} out
    } output {
        Array[File] out = glob("out*")
    } runtime {
        docker: "cowmoo/genomics-container:latest"
    }
}

task RSEMCalculateExpression {

    Array[File] reference_outs
    File fastq1
    File fastq2
    String sample_name

    command <<<
        outs=(${sep=" " reference_outs})
        for i in "$${"{"}outs[@]}"
        do
            bname=`basename $i`
            ln -s "$i" "$bname"
        done

        export PATH=/bowtie-1.3.1-linux-x86_64/:$PATH
        rsem-calculate-expression -p 30 --paired-end ${fastq1} ${fastq2} out ${sample_name}
    >>> output {
        File isoforms = sample_name + ".isoforms.results"
        File genes = sample_name + ".genes.results"
    } runtime {
        docker: "cowmoo/genomics-container:latest"
    }
}

task ExtractExons {

    File reference_gtf

    command {
        python3.9 <<CODE
        out = open('exon.gtf', "w")
        with open('${reference_gtf}') as f:
            while True:
                line = f.readline()

                if not line:
                    break

                if line.startswith('#'):
                    out.write(line)
                else:
                    fields = line.split("\t")
                    if fields[2] == "exon":
                        out.write(line)
        out.close()
        CODE
    } output {
       File out = "exon.gtf"
    } runtime {
        docker: "cowmoo/genomics-container:latest"
    }

}

workflow rnaseq {
    File fastq1
    File fastq2
    String sample_name

    File reference_gtf
    File reference_fa

    #call fastqc { input: fastq=fastq }

    call ExtractExons { input: reference_gtf=reference_gtf }

    call RSEMPrepareReference {
        input: reference_fa = reference_fa,
               exon_gtf = ExtractExons.out
    }

    call RSEMCalculateExpression {
        input: reference_outs = RSEMPrepareReference.out,
               fastq1 = fastq1,
               fastq2 = fastq2,
               sample_name = sample_name
    }

    #call STARAlignGenome {
    #    input: reference_fa = reference_fa,
    #           reference_gtf = reference_gtf
    #}

    #call STARAlignTranscripts {
    #    input: reference_outs = STARAlignGenome.out,
    #           fastq = fastq
    #}

    meta {
        author: "Paul Cao"
        email: "cowmoomoo@gmail.com"
        description: "## Bamstats \n This is the Bamstats workflow.\n\n Adding documentation improves clarity."
    }
}