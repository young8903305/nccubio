params.fqlist = "fqlist.txt"
params.fqpath = "$PWD/data/fq/"
params.grouplist = ""
params.genome = ""
params.gtf = "$PWD/data/genes.gtf"
params.fa = "$PWD/data/genome.fa"
params.p = 8
params.output = "results/"

/*
 * Initialization
 */

if( params.fqlist == "" ){
  error "--fqlist must be set."
}
else{
  fqlist_file = file(params.fqlist)
  if( !fqlist_file.exists() ){
    error "No fqlist file: ${fqlist_file}."
  }
  else{
    files_name = fqlist_file.text.split('\n|,')
    for( i = 0 ; i < files_name.length ; i = i + 1 ){
      current_file = file("${params.fqpath}${files_name[i]}")
      if( !current_file.exists() ){
        error "No fqfile file: ${current_file}."
      }
    }
  }
}

if( params.gtf == "" ){
  error "--gtf must be set."
}
else{
  gtf_file = file(params.gtf)
  if( !gtf_file.exists() ){
    error "No gtf file: ${gtf_file}."
  }
}

if( params.grouplist == "" ){
  error "--grouplist must be set."
}
else{
  grouplist_file = file(params.grouplist)
  if( !grouplist_file.exists() ){
    error "No grouplist file: ${grouplist_file}."
  }
}

if( params.fa == "" ){
  error "--fa must be set."
}
else{
  fa_file = file(params.fa)
  if( !fa_file.exists() ){
    error "No fa file: ${fa_file}."
  }
}

if( params.genome == "" ){
  error "--genome must be set."
}
else{
  genome_folder = file(params.genome)
  if( !genome_folder.exists() ){
    error "No genome folder: ${genome_folder}."
  }
}

if( params.output == "" ){
  error "--output must be set."
}
else{
  output_folder = file(params.output)
  if( !output_folder.exists() ){
    output_folder.mkdirs()
  }
  assemblies_file = file("${output_folder}/assemblies.txt")
}

Channel
  .fromPath(fqlist_file)
  .splitCsv(header:["fq1","fq2"])
  .set { fq_pairs }

/*
 * Align the RNA-seq reads to the genome
 * 1| Map the reads for each sample to the reference genome.
 * Assemble expressed genes and transcripts
 * 2| Assemble transcripts for each sample.
 */

process AlignAndAssemble {

  maxForks 1

  input:
    set ( fq1, fq2 ) from fq_pairs

  output:
    set ( state, id, file('thout_folder'), file('clout_folder') ) into thcls
    file thout_folder
    file clout_folder

  script:
    id = fq1.split('_')[0]
    state = ""
    if( file("${output_folder}/${id}_thout").exists() & file("${output_folder}/${id}_clout").exists() ){
      state = "exists"
    }
    else{
      state = "execute"
      fq1 = file("${params.fqpath}${fq1}")
      fq2 = file("${params.fqpath}${fq2}")
    }

    if( state == "exists" ){
      """
      echo "exists" > thout_folder
      echo "exists" > clout_folder
      """
    }
    else{
      """
      tophat -p $params.p -G $gtf_file -o thout_folder $genome_folder/genome $fq1 $fq2
      cufflinks -p $params.p -o clout_folder thout_folder/accepted_hits.bam
      """
    }
}

/*
 * 3| Create a file called assemblies.txt that lists the assembly file for each sample.
 */
process CreateAssemblyFile{
  input:
    set ( state, id, thout_folder, clout_folder ) from thcls

  output:
    file assembly_file

  script:
    if( state == "execute" ){
      thout_folder.moveTo("${output_folder}/${id}_thout")
      clout_folder.moveTo("${output_folder}/${id}_clout")
    }

    """
    touch assemblies_file
    echo ${output_folder}/${id}_clout/transcripts.gtf > assembly_file
    """
}

assemblies_file = assembly_file.collectFile( name:"${output_folder}/assemblies.txt")

/*
 * 4| Run Cuffmerge on all your assemblies to create a single merged transcriptome annotation.
 */

process Cuffmerge {

  input:
    file assemblies_file

  output:
    set ( state, file(merged_asm) ) into merge_res

  script:
    state = ""
    if( file("${output_folder}/merged_asm").exists() ){
      state = "exists"
    }
    else{
      state = "execute"
    }

    if( state == "exists" )
      """
      echo "exists" > merged_asm
      """
    else
      """
      cuffmerge -g $gtf_file -s $fa_file -p $params.p $assemblies_file
      """

}

/*
 * Identify differentially expressed genes and transcripts
 * 5| Run Cuffdiff by using the merged transcriptome assembly along with the BAM files from TopHat for each replicate.
 */

process Cuffdiff {

  input:
    set ( state, merged_asm ) from merge_res
  output:
    file diff_out

  script:
    if( state == "execute" ){
      merged_asm.moveTo("${output_folder}/merged_asm")
    }
    merged_gtf = file("${output_folder}/merged_asm/merged.gtf")
    groups = grouplist_file.text.split('\n')
    group_str = ""
    for( i=0 ; i<groups.length ; i=i+1 ){
      entries = groups[i].split(',')
      for( j=0 ; j<entries.length ; j=j+1 ){
        group_str = "${group_str}${output_folder}/${entries[j]}_thout/accepted_hits.bam"
        if( j!=entries.length-1 ){
          group_str = "${group_str},"
        }
      }
      group_str = "${group_str} "
    }
    """
    cuffdiff -o diff_out -b $fa_file -p $params.p -u $merged_gtf ${group_str}
    """

}

diff_out.subscribe{
  it.moveTo("${output_folder}/diff_out")
}
