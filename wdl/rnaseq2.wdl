task SalmonQuantify {
    File fastq1
    File fastq2
    String input_sample
    String input_condition
    Array[File] index

    command <<<
        mkdir index
        outs=(${sep=" " index})
        for i in "$${"{"}outs[@]}"
        do
            bname=`basename $i`
            ln -s "$i" index/"$bname"
        done

        /salmon-1.9.0_linux_x86_64/bin/salmon quant -i index -l A \
         -1 ${fastq1} \
         -2 ${fastq2} \
         -p 8 --validateMappings -o quant

        mv quant/quant.sf ${input_sample}
    >>>
    output {
        File quant = "${input_sample}"
        String sample = input_sample
        String condition = input_condition
    }
    runtime {
        docker: "cowmoo/genomics-container:latest"
    }
}

task SalmonIndex {
    File transcriptome

    command {
        /salmon-1.9.0_linux_x86_64/bin/salmon index -t ${transcriptome} -i transcript_index
    }
    output {
        Array[File] index = glob("transcript_index/*.*")
    } runtime {
        docker: "cowmoo/genomics-container:latest"
    }
}

task DESeq2Report {
    Array[File] quants
    File sample_sheet

    command <<<
        outs=(${sep=" " quants})
        for i in "$${"{"}outs[@]}"
        do
            bname=`basename $i`
            mkdir /"$bname"
            ln -s "$i" /"$bname"/quant.sf
        done

        sed -i '1s/^/Sample\tfq1\tfq2\tCondition\n/' ${sample_sheet}

        ln -s ${sample_sheet} /samples.txt

        Rscript -e "rmarkdown::render('/genomics-portfolio/reports/sample_report.Rmd', output_dir='.')"
    >>> output {
        File pdf = "sample_report.pdf"
    } runtime {
        docker: "cowmoo/genomics-container:latest"
    }
}

workflow rnaseq {
    File sample_sheet
    File transcriptome

    Array[Array[String]] samples = read_tsv(sample_sheet)

    call SalmonIndex { input: transcriptome=transcriptome }

    scatter (sample in samples) {
        call SalmonQuantify {
            input: index = SalmonIndex.index,
                   fastq1 = sample[1],
                   fastq2 = sample[2],
                   input_sample = sample[0],
                   input_condition = sample[3]
        }
    }

    call DESeq2Report {
        input: sample_sheet = sample_sheet,
               quants = SalmonQuantify.quant
    }

    meta {
        author: "Paul Cao"
        email: "cowmoomoo@gmail.com"
        description: "## Bamstats \n This is the Bamstats workflow.\n\n Adding documentation improves clarity."
    }
}