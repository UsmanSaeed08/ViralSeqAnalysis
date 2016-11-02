package positiveInsideRule;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import org.biojava.bio.BioException;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojavax.SimpleNamespace;
import org.biojavax.bio.seq.RichSequence;
import org.biojavax.bio.seq.RichSequenceIterator;


public class Fasta_Parser2 {

	private ArrayList <Elements_file> parsed_file;
	private String file_path;
	private String temp_file[][];
	private String temp_file_final[][];
	private int temp_file_index;
	public static final HitStats_holder holder = new HitStats_holder();
	public static long numberOfprots = 0 ;

	
	public Fasta_Parser2(ArrayList <Elements_file> p, String path){
		parsed_file = p;
		file_path = path;
		temp_file = new String[1000][2];
		temp_file_final = new String[1000][2];
		this.temp_file_index = 0;
		
	}

	public ArrayList <Elements_file> parser(){	//reads fasta file and adds in holder

		BufferedReader br ;
		String acc = "";
		
		//ArrayList <Elements_file> temp_ = new ArrayList <Elements_file>(); 

		try{
			
			//System.out.print(file_path+"\n");
			
			br = new BufferedReader(new FileReader(file_path)); //add file path here
		    
		    Alphabet alpha = AlphabetManager.alphabetForName("PROTEIN");
		    SimpleNamespace ns = new SimpleNamespace("biojava");
		    RichSequenceIterator iterator = RichSequence.IOTools.readFasta(br,
		    alpha.getTokenization("token"), ns);
		    
		    
		    int i = 0;
		    while (iterator.hasNext()) {	//populate temp_file
		        RichSequence rec = iterator.nextRichSequence();
		        acc = rec.getAccession();
		        String seq = rec.seqString();
		        
		        	temp_file[i][0] = acc;
		        	temp_file[i][1] = seq;
		        
		        i++;       
		    }
		    int z = 0;	//index for temp_file_final
		    int j =0;
		    while(j<i){	//while end of temp_file
		    	if(!temp_file[j][0].contains("orf") && temp_file[j][0].contains("-")){	//if it is not an orf sequence then----for each entry in file--make elements like before but with 
		    											// fused sequences which is done through loops
		    		int no = 0;
		    		
	    			int idx = temp_file[j][0].indexOf("-");
	    			char n = temp_file[j][0].charAt(idx+1);
	    			no = Character.getNumericValue(n);
	    			int k =1;
	    			//int x = j;	//index to use for temp_file
	    			String ac_cut = temp_file[j][0].substring(0, idx);
	    			String ac = make_ac(idx,temp_file[j][0]);	//this fuction makes new acc similar to one in phobius results
	    			String s = makeString(j,ac_cut);
		    		temp_file_final[z][0] = ac;
			        temp_file_final[z][1] = s;
			        
			        j = this.temp_file_index;
			        z++;
	    			/*
		    			while(k<=no){
		    				if(k==1){
		    					temp_file_final[z][0] = ac;
			    				temp_file_final[z][1] = temp_file[x][1];
		    				}
		    				else if(k>1){
			    				temp_file_final[z][1].concat(temp_file[x][1]);
		    				}
		    				
		    				k++;
		    				x++;
		    			}
		    			//System.out.print(no);
		    			j = x;
		    			z++;
		    			*/
		    		
		    	}
		        else{
		        	//String temp = acc.s
		        	temp_file_final[z][0] = temp_file[j][0];
			        temp_file_final[z][1] = temp_file[j][1];
			        j++;
			        z++;
		        	
		        }
		    }
		    
		    //Now put it in element file to return
		    Elements_file element = new Elements_file();
		    //int t = 0;
		    for(int x=0; x<=z-1;x++){		//x same because sequences are parallel in parsed file and temp file final
		    	element = search_parse(temp_file_final[x][0]);
		    	 if (element!=null && x<parsed_file.size() && element.tm_no>0){
				        if (element.file_name == parsed_file.get(x).file_name && element.id.contains(temp_file_final[x][0])){
				        	
				        	//System.out.print(element.id + " --> To Do\n");
				        	element.seq = temp_file_final[x][1];
				        	element = get_tm_section(element);
				        	numberOfprots++;
				        	parsed_file.set(x, element);		        	
				        }
				    }
		    }
		}
		catch(IOException e){
			e.printStackTrace();
			System.exit(0);
		} catch (BioException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			//System.out.print("error at " + acc);
			System.exit(0);
		}
		return parsed_file;
	}
	
	private String makeString(int x, String ac) {
		String seq = "";
		//System.out.print(ac+"\n");
		while(temp_file[x][0].contains(ac)){
			seq = seq + temp_file[x][1];
			//seq.concat(temp_file[x][1]);
			x++;
			if( temp_file[x+1][0] == null){
				break;
			}
		}
		
		temp_file_index = x;

		return seq;
	}

	private String make_ac(int idx, String accession) {
		String temp1 = accession.substring(idx+2);
		String temp2 = accession.substring(0, idx);
		return temp2 + temp1;
	}

	private Elements_file search_parse(String acc) {
		for (int i =0; i<=parsed_file.size()-1; i++){
			if (parsed_file.get(i).id.contains(acc)){
				//System.out.print(acc+"\n");
				return parsed_file.get(i);
			}
		}
		// TODO Auto-generated method stub
		return null;
	}

	private Elements_file get_tm_section(Elements_file e){
		// 
		for (int i =0; i<= e.tm_no-1; i++){
			String tm = e.seq.substring(e.start[i]-1, e.end[i]+1);	//-2 and +2 to get the flank
			e.tm_segments.add(tm);
		}
		if(!e.seq.isEmpty()){
			getSegmentsAndStats(e);
		}
		
		return e;
		
	}

	private void getSegmentsAndStats(Elements_file e) {
		// TODO Auto-generated method stub
		// get all the segments
		ArrayList<String> Segments = new ArrayList<String>();
		String tempSeq = e.seq;
		for(int j =0; j<=e.tm_no-1;j++){
			int st = e.start[j]; 
			int end = e.end[j];

			if(j==0){
				// get the first part
				String part = tempSeq.substring(0, st-1);
				Segments.add(part);
				// get the tm part
				part = tempSeq.substring(st, end);
				Segments.add(part);
				if (j+1 >e.tm_no-1){
					// for the last segment
					int len = tempSeq.length();
					part = tempSeq.substring(end, len);
					Segments.add(part);
				}
			}
			else{
				//int stLast = e.start[j-1]; 
				int endLast = e.end[j-1];
				// get the first part
				String part = tempSeq.substring(endLast+1, st-1);
				Segments.add(part);
				// get the tm part
				part = tempSeq.substring(st, end);
				Segments.add(part);
				if (j+1 >e.tm_no-1){
					// for the last segment
					int len = tempSeq.length();
					part = tempSeq.substring(end, len);
					Segments.add(part);
				}
			}

		}
		// since now the tm segment positions are known and also the topolgy
		// should now calculate the statistics

		ArrayList<String> InsideSegments = new ArrayList<String>();
		ArrayList<String> OutsideSegments = new ArrayList<String>();
		ArrayList<String> TmSegments = new ArrayList<String>();
		
		
		if (e.nterm == 1){	// first is Inside
			int lastEntry = 0; // 1 for inside, 2 for TM and 3 for outside
			for (int j = 0; j<= Segments.size()-1;){
				InsideSegments.add(Segments.get(j));
				j = j+1;
				lastEntry = 1 ; // For inside identifier
				if(j<=Segments.size()-1){
					TmSegments.add(Segments.get(j));
					j = j+1;
					lastEntry = 2 ; // For TM identifier
				}
				if(j<=Segments.size()-1){
					OutsideSegments.add(Segments.get(j));
					j = j+1;
					lastEntry = 3 ; // For outside identifier
				}
				if(j<=Segments.size()-1){
					TmSegments.add(Segments.get(j));
					j = j+1;
					lastEntry = 2 ; // For TM identifier
				}
			}
			// e.cterm 1 for inside and 0 for outside
			if (e.cterm == 1 && lastEntry != 1){
				// bad
				// no use of the whole protein here.. so it is useless to count all the hit stats, so just continue
				return;
			}
			else if(e.cterm == 0 && lastEntry != 3){
				// bad
				// no use of the whole protein here.. so it is useless to count all the hit stats, so just continue
				return;
			}
		}
		else if (e.nterm == 0){				// first is Outside
			int lastEntry = 0; // 1 for inside, 2 for TM and 3 for outside
			for (int j = 0; j<= Segments.size()-1;){
				OutsideSegments.add(Segments.get(j));
				j = j+1;
				lastEntry = 3;
				if(j<=Segments.size()-1){
					TmSegments.add(Segments.get(j));
					j = j+1;
					lastEntry = 2;
				}

				if(j<=Segments.size()-1){
					InsideSegments.add(Segments.get(j));
					j = j+1;
					lastEntry = 1;
				}

				if(j<=Segments.size()-1){
					TmSegments.add(Segments.get(j));
					j = j+1;
					lastEntry = 2;
				}
			}
			if (e.cterm == 1 && lastEntry != 1){
				// bad
				// no use of the whole protein here.. so it is useless to count all the hit stats, so just continue
				return;
			}
			else if(e.cterm == 0 && lastEntry != 3){
				// bad
				// no use of the whole protein here.. so it is useless to count all the hit stats, so just continue
				return;
			}
		}
		Segments = null;
		HitStats_holder thisprotein = new HitStats_holder();
		
		for (int x = 0;x <= InsideSegments.size()-1; x++){
			thisprotein.obj_hitstats_inside = Check_hit_tm.count(InsideSegments.get(x), thisprotein.obj_hitstats_inside);
			thisprotein.obj_hitstats_inside.NumberOfthisSegment = thisprotein.obj_hitstats_inside.NumberOfthisSegment + 1 ;
		}
		for (int x = 0;x <= TmSegments.size()-1; x++){
			// commented out below inorder to take one periplasm and one cytoplasm
			thisprotein.obj_hitstats = Check_hit_tm.count(TmSegments.get(x), thisprotein.obj_hitstats);
			thisprotein.obj_hitstats.NumberOfthisSegment = thisprotein.obj_hitstats.NumberOfthisSegment + 1 ;
			
		}
		for (int x = 0;x <= OutsideSegments.size()-1; x++){
			thisprotein.obj_hitstats_outsdide = Check_hit_tm.count(OutsideSegments.get(x), thisprotein.obj_hitstats_outsdide);
			thisprotein.obj_hitstats_outsdide.NumberOfthisSegment = thisprotein.obj_hitstats_outsdide.NumberOfthisSegment + 1 ;
		}
				
		// now the array has percentages and indivisual integers have count
		// now add thisprotein to holder because otherwise all the percentages and counts would
		// get mixed up
		holder.obj_hitstats_inside = Check_hit_tm.SumUp(holder.obj_hitstats_inside,thisprotein.obj_hitstats_inside);
		holder.obj_hitstats = Check_hit_tm.SumUp(holder.obj_hitstats,thisprotein.obj_hitstats);
		holder.obj_hitstats_outsdide = Check_hit_tm.SumUp(holder.obj_hitstats_outsdide,thisprotein.obj_hitstats_outsdide);

	}

}
