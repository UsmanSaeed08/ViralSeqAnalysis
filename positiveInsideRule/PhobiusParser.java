package positiveInsideRule;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

public class PhobiusParser {

	private String path;
	private int st_tm;
	private int end_tm;
	private ArrayList <Elements_file> parsed_file;
	String name; //file name

	public PhobiusParser(String p, String n){
		path =p;
		st_tm = 0;
		end_tm = 0;
		name = n;
	}
	
	public ArrayList<Elements_file> open_result(){

		BufferedReader br = null;
		BufferedReader br2 =null;

		String sCurrentLine;
		String subLine = "";
		String lastLine = "";
		String cross = "//";
		parsed_file = new ArrayList<Elements_file>();

		int line_number = 0;

		try{ 
			br = new BufferedReader(new FileReader(path));

			while ((sCurrentLine = br.readLine()) != null) {//reads whole file
				
				String id =null;
				if (sCurrentLine.charAt(0) == 'I' && sCurrentLine.charAt(1) == 'D'){
					//get id
					br2 = new BufferedReader(new FileReader(path));
					Elements_file element = new Elements_file();
					//id = get_id(sCurrentLine);
					id = sCurrentLine;
					element.id = id;
					element.file_name = name;
					int index = 0;
				

					for (int x = 0; x<=line_number; x++){
						subLine = br2.readLine();		//to reach the current line, current line is the id line
					}
					lastLine = subLine;

					do{
						subLine = br2.readLine();
						if (subLine!= null && subLine.length()>5){
				
							//if (subLine.charAt(5) == 'T'){	//change to S for singal pep
							/*
							if (subLine.charAt(5) == 'S'){	//change to S for singal pep
								//get tm segment
								SP_segment(subLine);
								if(element.SP == true){
									System.out.print("\n ERROR--- SP already exists \n");
									System.exit(0);
								}
								element.SP = true;
								element.sp_begin = st_tm;
								element.sp_end = end_tm;
								
							}*/
							//else if (subLine.charAt(5) == 'T'){	//change to S for singal pep
							if (subLine.charAt(5) == 'T'){	//change to S for singal pep
								if (element.tm_no == 0){	// first TM found, so can get the starting terminal of protein
									// i need the last line processed over here
									// got the lastLine
									// if line does not contain NON and contains CYTOPLASMIC. then it is inside
									// else if line contains NON CYTOPLASMIC. then it is outside
									if (!(lastLine.contains("NON")) && (lastLine.contains("CYTOPLASMIC."))){
										element.nterm = 1;
									}
									else if(lastLine.contains("NON CYTOPLASMIC.")){
										element.nterm = 0;
									}
								}
								//get tm segment
								tm_segment(subLine);
								element.tm_no = element.tm_no + 1;
								element.start[index] = st_tm;
								element.end[index] = end_tm;
								index ++;
							}
						}
						if(subLine != null && subLine.charAt(0) != 'I' && subLine.charAt(1) != 'D' && !subLine.contains("//")){
							lastLine = subLine;
						}
					}
					while(subLine != null && subLine.charAt(0) != 'I' && subLine.charAt(1) != 'D');
					br2.close();
					// So when it comes out of the loop. the last line has not been changed to the ID line
					if (!(lastLine.contains("NON")) && (lastLine.contains("CYTOPLASMIC."))){
						element.cterm = 1;
					}
					else if(lastLine.contains("NON CYTOPLASMIC.")){
						element.cterm = 0;
					}
					parsed_file.add(element);
				}
				line_number ++;
			}
			br.close();
			
			/*
			  	for (int i =0; i<=parsed_file.size()-1; i++){
				System.out.print("ID: " + parsed_file.get(i).id + "\n");
				System.out.print("TM_NO: " + parsed_file.get(i).tm_no + "\n");
				for (int j = 0; j<parsed_file.get(i).start.length-1;j++){
					System.out.print("Start: " + parsed_file.get(i).start[j]);
					System.out.print("    End: " + parsed_file.get(i).end[j]);
					System.out.print("\n");
				}
				System.out.print("\n\n\n");
			}*/
		}
		
		catch(IOException e){
			e.printStackTrace();
		}
		return parsed_file;
	}

	private void tm_segment (String st1){

		st1 = st1.substring(13, st1.length());
		st1 = st1.trim();
		int mid = st1.length()/2;

		String st2 = st1.substring(mid);
		st1 = st1.substring(0, mid);
		st2 = st2.trim();
		st1 = st1.trim();
		st_tm = Integer.valueOf(st1);
		end_tm = Integer.valueOf(st2);

	}
	private void SP_segment (String st1){

		st1 = st1.substring(12, st1.length());
		st1 = st1.trim();
		int mid = st1.length()/2;

		String st2 = st1.substring(mid);
		st1 = st1.substring(0, mid);
		st2 = st2.trim();
		st1 = st1.trim();
		st_tm = Integer.valueOf(st1);
		end_tm = Integer.valueOf(st2);

	}

	private String get_id (String line){
		int start =0;
		int end = 0;
		for (int i =0; i<line.length()-1; i++){
			if (line.charAt(i) == 'g'){
				start=i;
			}
			if (line.charAt(i) == '|' && i>18){
				end =i;
			}
		}
		String id = line.substring(start, end);
		return id;
	}



}
