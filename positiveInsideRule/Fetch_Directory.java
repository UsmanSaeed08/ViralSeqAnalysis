package positiveInsideRule;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.data.general.DefaultPieDataset;

import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.category.DefaultCategoryDataset;

public class Fetch_Directory {

	String path;
	int count; 

	private ArrayList<String> files_results; //List of files in directory...use it to send to parser which would parse each file...
	// or loop here to call parser every time
	//public Elements_File_Holder holder[];


	public Fetch_Directory(){

	}

	public Fetch_Directory(String p){
		path = p;
		final File folder = new File(path);
		listFilesForFolder(folder);
		count = 0;
	}



	public void listFilesForFolder(final File folder) {

		int count = 0; 
		files_results = new ArrayList<String>();

		for (final File fileEntry : folder.listFiles()) {
			if (fileEntry.isDirectory()) {
				listFilesForFolder(fileEntry);
			} else {
				count++;
				files_results.add(fileEntry.getName().toString());
			}
		}

		//for (int index=0; index<=files_results.size()-1;index++){	
		//System.out.print(files_results.get(index)+"\n");
		//}
		System.out.println("\n total file " + count + " \n");
		// all file names received
		// call files_handler
		files_handler();
	}


	public void files_handler(){

		String name_result_file;
		//ArrayList <Elements_file> temp_ = new ArrayList <Elements_file>(); 
		Elements_File_Holder holder[] = new Elements_File_Holder[1030];

		String pfasta = "F:/virus_project/input/genomes.domains";
		//String pfasta = "F:/virus_project/1/faa";

		for (int index=0; index<=files_results.size()-1;index++){	
			name_result_file = files_results.get(index);
			// call phobius parser			
			PhobiusParser obj = new PhobiusParser(path + "/" + name_result_file, name_result_file);

			holder[index] = new Elements_File_Holder();						//Extracting the coordinates of TM
			holder[index].file_name = name_result_file;
			holder[index].obj = obj.open_result();

		}
		int index =0;
		while(index<=files_results.size()-1){
			name_result_file = files_results.get(index);
			name_result_file = name_result_file.replaceAll(".tm", ".faa");
			//Extracting the TM through coordinates
			Fasta_Parser2 obj1 = new Fasta_Parser2(holder[index].obj,pfasta + "/" + name_result_file);
			holder[index].obj = obj1.parser();
			//holder[index].obj = getscore((holder[index].obj)); // checking the miss lead in tm extracted and tm count
			
			index++;
		}

		System.out.print("\nDone -> "+ count+"\n" );
		System.out.print("Done -> seq extraction and now BLAST \n");
		printResults(Fasta_Parser2.holder, Fasta_Parser2.numberOfprots);
		/*
		index =0;

		//HitStats obj_hitstats = new HitStats();	//for tm section only
		HitStats_holder rule = new HitStats_holder();
		
		int queryseqcount =0;
		int totalseq =0;

		while(index<=files_results.size()-1){	//go through all the files
			ArrayList <Elements_file> p = holder[index].obj;
			index++;
			Check_hit_tm oj = new Check_hit_tm();


			for (int i =0; i<= p.size()-1;i++){
				totalseq++;
				
				//if (p.get(i).tm_no>0 && p.get(i).SP == true){
				if (p.get(i).tm_no>0){
					queryseqcount++;
					//obj_hitstats = oj.getcounts(p.get(i).seq,obj_hitstats,p.get(i).start,p.get(i).end,p.get(i).tm_no);//Commenting out cuz dont need aa distribution
					//obj_hitstats = oj.getcounts_sp(p.get(i).seq,obj_hitstats,p.get(i).sp_begin,p.get(i).sp_end);
					if (p.get(i).tm_no == 1){
						rule = oj.getcounts_rule_one(p.get(i).seq,rule,p.get(i).start,p.get(i).end,p.get(i).tm_no,p.get(i));
					}
					else{
						rule = oj.getcounts_rule(p.get(i).seq,rule,p.get(i).start,p.get(i).end,p.get(i).tm_no);
					}
					rule.obj_hitstats.hits_no++;	
				}
			}
		}
//TM
		System.out.print("\n TM \n");
		for(int i =0; i<= rule.obj_hitstats.aa_seq.length()-1;i++){
			System.out.print("\n"+ rule.obj_hitstats.aa_seq.charAt(i) + " --- " + rule.obj_hitstats.aa[i]);
		}
		System.out.print("\n INSIDE \n");
//Inside		
		for(int i =0; i<= rule.obj_hitstats.aa_seq.length()-1;i++){
			System.out.print("\n"+ rule.obj_hitstats.aa_seq.charAt(i) + " --- " + rule.obj_hitstats_inside.aa[i]);
		}
		System.out.print("\n OUTSIDE \n");
//Outside
		for(int i =0; i<= rule.obj_hitstats.aa_seq.length()-1;i++){
			System.out.print("\n"+ rule.obj_hitstats.aa_seq.charAt(i) + " --- " + rule.obj_hitstats_inside.aa[i]);
		}
		System.out.print("\n \n");
		
		System.out.print("\n"+"no of TM: " + queryseqcount);
		System.out.print("\n"+"no of Total seq: " + totalseq);

*/



		/*
		HitStats obj_hitstats = new HitStats();	//origin of hitstats obj...cuz here is where all sequences for all files would be done for BLAST

		//*************** IMPORTANT**********************
		int SequencesCount = 0;	//number of total query sequences with TM and would be used for BLAST
								//number of total hits are in HitStats 

		//HitStats hitstatObj = new HitStats();	//for the Query case

		while(index<=files_results.size()-1){	// while all files--->
			/*	************* THE NON VIRAL ALL HOMOLOG *********
			Simap_client blast_obj = new Simap_client(holder[index].obj,obj_hitstats); // for each file run blast
			holder[index].obj = blast_obj.run();
			obj_hitstats = blast_obj.hitstatObj;
			SequencesCount = SequencesCount + blast_obj.queryseqcount;
			index++;
			// ******** THE QUERY CASE *************

			// for each seq in this file if tm seq then call and calc

			ArrayList <Elements_file> p = holder[index].obj;
			Check_hit_tm oj = new Check_hit_tm();

			int queryseqcount =0;
			for (int i =0; i<= p.size()-1;i++){
				if (p.get(i).tm_no>0){
					queryseqcount++;
					obj_hitstats = oj.getcounts(p.get(i).seq,obj_hitstats,p.get(i).start,p.get(i).end,p.get(i).tm_no);//was in if above
					obj_hitstats.hits_no++;	
			}

		}
			SequencesCount = SequencesCount + queryseqcount;
			index++;
		}
			//************** QUERY CASE ENDS*************
		makechart(obj_hitstats);
		System.out.print("\n"+"no of TM Query: " + SequencesCount);


		System.out.print("\nDone -> BLAST \n");
		 */
	}
	
	private void printResults(HitStats_holder holder, long n){
		System.out.print("\n Printing Results for TM: "+"\n");
		for(int i =0; i<= holder.obj_hitstats.aa_seq.length()-1;i++){
			//System.out.print("\n"+ holder.obj_hitstats.aa_seq.charAt(i) + " --- " + holder.obj_hitstats.aa[i]);
			//System.out.print(holder.obj_hitstats.aa[i]/numberOfprots+"\n");
			System.out.print(holder.obj_hitstats.aa[i]/holder.obj_hitstats.totallength+"\n");
		}
		
		
		
		System.out.print("\n Printing Results for Inside: "+"\n");
		for(int i =0; i<= holder.obj_hitstats_inside.aa_seq.length()-1;i++){
			//System.out.print("\n"+ holder.obj_hitstats_inside.aa_seq.charAt(i) + " --- " + holder.obj_hitstats_inside.aa[i]);
			//System.out.print(holder.obj_hitstats_inside.aa[i]/numberOfprots+"\n");
			System.out.print(holder.obj_hitstats_inside.aa[i]/holder.obj_hitstats_inside.totallength+"\n");
		}
		
		
		System.out.print("\n Printing Results for Outside: "+"\n");
		for(int i =0; i<= holder.obj_hitstats_outsdide.aa_seq.length()-1;i++){
			//System.out.print("\n"+ holder.obj_hitstats_outsdide.aa_seq.charAt(i) + " --- " + holder.obj_hitstats_outsdide.aa[i]);
			//System.out.print(holder.obj_hitstats_outsdide.aa[i]/numberOfprots+"\n");
			System.out.print(holder.obj_hitstats_outsdide.aa[i]/holder.obj_hitstats_outsdide.totallength+"\n");
		}
		/*
		System.out.print("hydrophobic in tm: "+holder.obj_hitstats.hydrophobic/numberOfprots+"\n");
		System.out.print("hydrophilic in tm: "+holder.obj_hitstats.hydrophilic/numberOfprots+"\n");
		System.out.print("KnR in tm: "+holder.obj_hitstats.KnR/numberOfprots+"\n");
		
		System.out.print("hydrophobic Inside: "+holder.obj_hitstats_inside.hydrophobic/numberOfprots+"\n");
		System.out.print("hydrophilic Inside: "+holder.obj_hitstats_inside.hydrophilic/numberOfprots+"\n");
		System.out.print("KnR Inside: "+holder.obj_hitstats_inside.KnR/numberOfprots+"\n");
		
		System.out.print("hydrophobic Outside: "+holder.obj_hitstats_outsdide.hydrophobic/numberOfprots+"\n");
		System.out.print("hydrophilic Outside: "+holder.obj_hitstats_outsdide.hydrophilic/numberOfprots+"\n");
		System.out.print("KnR Outside: "+holder.obj_hitstats_outsdide.KnR/numberOfprots+"\n");
		*/
		System.out.print("hydrophobic in tm: "+holder.obj_hitstats.hydrophobic/holder.obj_hitstats.totallength+"\n");
		System.out.print("hydrophilic in tm: "+holder.obj_hitstats.hydrophilic/holder.obj_hitstats.totallength+"\n");
		System.out.print("KnR in tm: "+holder.obj_hitstats.KnR/holder.obj_hitstats.totallength+"\n");
		
		System.out.print("hydrophobic Inside: "+holder.obj_hitstats_inside.hydrophobic/holder.obj_hitstats_inside.totallength+"\n");
		System.out.print("hydrophilic Inside: "+holder.obj_hitstats_inside.hydrophilic/holder.obj_hitstats_inside.totallength+"\n");
		System.out.print("KnR Inside: "+holder.obj_hitstats_inside.KnR/holder.obj_hitstats_inside.totallength+"\n");
		
		System.out.print("Positive Residues Inside: "+holder.obj_hitstats_inside.KnR+"\n");
		System.out.print("Total Residues Inside: " + holder.obj_hitstats_inside.totallength+"\n");
		
		System.out.print("hydrophobic Outside: "+holder.obj_hitstats_outsdide.hydrophobic/holder.obj_hitstats_outsdide.totallength+"\n");
		System.out.print("hydrophilic Outside: "+holder.obj_hitstats_outsdide.hydrophilic/holder.obj_hitstats_outsdide.totallength+"\n");
		System.out.print("KnR Outside: "+holder.obj_hitstats_outsdide.KnR/holder.obj_hitstats_outsdide.totallength+"\n");
		
		System.out.print("Positive Residues Outside: "+holder.obj_hitstats_outsdide.KnR+"\n");
		System.out.print("Total Residues Outside: "+ holder.obj_hitstats_outsdide.totallength+"\n");
		

		System.out.print("No of Prots "+n+"\n");

	}

	public void makechart(HitStats hits){
		for(int i =0; i<= hits.aa_seq.length()-1;i++){
			System.out.print("\n"+ hits.aa_seq.charAt(i) + " --- " + hits.aa[i]);
		}
		System.out.print("\n"+"no of Query: " + hits.hits_no);


		String dir = "D:/virus_project";
		DefaultPieDataset pieDataset = new DefaultPieDataset();
		DefaultCategoryDataset dataset = new DefaultCategoryDataset();

		for(int i =0; i<=hits.aa_seq.length()-1;i++){
			pieDataset.setValue(Character.toString(hits.aa_seq.charAt(i)), new Float( hits.aa[i]));
			dataset.setValue(hits.aa[i], "No", Character.toString(hits.aa_seq.charAt(i)));
		}

		JFreeChart chart = ChartFactory.createPieChart("AminoAcid Distribution",pieDataset,true,true,false);
		JFreeChart chart_bar = ChartFactory.createBarChart("AminoAcid Distribution Bargraph","aa", "No", dataset, PlotOrientation.VERTICAL, false,true, false);

		try {

			//ChartUtilities.saveChartAsJPEG(new File(dir+"\\"+ "hits_" +"piechart_QueryStat.jpg"), chart, 500,300);
			ChartUtilities.saveChartAsJPEG(new File(dir+"\\"+ "piechart_EColiTmStat.jpg"), chart, 500,300);
			ChartUtilities.saveChartAsJPEG(new File(dir+"\\"+ "bargraph_EColiTmStat.jpg"), chart_bar, 500,300);

		} catch (Exception e) {

			System.out.println("Problem occurred creating chart.");

		}
	}

	private ArrayList <Elements_file> getscore(ArrayList <Elements_file> p) {
		// TODO Auto-generated method stub
		int index = 0;
		//float scorep1, scorep2;
		//float scoren1, scoren2;
		int flag =0;
		while(index < p.size()-1){	// loop to go through all the sequences in one file

			if (p.get(index).tm_no != p.get(index).tm_segments.size()){
				if (flag == 0){
					System.out.print("File name: "+p.get(index).file_name + "\n");
					flag =1;
				}

				System.out.print(p.get(index).id + "\n");
				count ++;
			}index++;

			/*
			if(p.get(index).tm_no > 0){
				if (p.get(index).tm_no == p.get(index).tm_segments.size()){
				int v= p.get(index).tm_no;
				p.get(index).intializeScore(v);
				for (int i = 0; i<=v-1;i++){	//loop to go through all the tm segments in that sequence
					scorep1 = 0;
					scoren1 = 0;
					String s = p.get(index).tm_segments.get(i).toString();	
					for (int j =0; j<= s.length()-1; j++){		//loop to see each character in tm segment
						if (s.charAt(j)=='K' || s.charAt(j)=='R' || s.charAt(j)=='H'){
							scorep1 = scorep1 + 1;
						}

						if (s.charAt(j)=='D' || s.charAt(j)=='E'){
							scoren1 = scoren1 + 1;
						}
					}
					scorep1 = scorep1/s.length();	//score for segment calculated
					scoren1 = scoren1/s.length();
					p.get(index).segment_scorep[i] = scorep1;
					p.get(index).segment_scoren[i] = scoren1;
				}
			}
		}
			index ++; */
		}


		return p;
	}



	private int[] tm_exists(Elements_File_Holder holder[]){
		int i =0;
		int file_to_Check[] = new int[1000];
		for (i=0;i<999;i++)
			file_to_Check[i] =0;
		i = 0;
		int j = 0;
		int z = 0;
		while (holder[i]!= null){
			while (!holder[i].obj.isEmpty()){
				if (holder[i].obj.get(j).tm_no>0){
					file_to_Check[z] = i;
					z++;
					break;
				}
				j++;
			}
			i ++;
			j =0;
		}
		return file_to_Check;
	}



}
