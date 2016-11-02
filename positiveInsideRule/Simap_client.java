package positiveInsideRule;

import java.util.ArrayList;

import org.apache.commons.lang3.StringUtils;

import mips.gsf.de.simapclient.client.SimapAccessWebService;
import mips.gsf.de.simapclient.datatypes.HitSet;
import mips.gsf.de.simapclient.datatypes.ResultParser;



// simap client for general run

public class Simap_client {

	private ArrayList <Elements_file> parsed_file;
	private int rangeIndex;
	public HitStats hitstatObj;
	public int queryseqcount;

	public Simap_client(ArrayList <Elements_file> p,HitStats h){
		parsed_file = p;
		hitstatObj = h;
		rangeIndex = 0;
		queryseqcount = 0;

	}

	public  ArrayList <Elements_file> run(){
		//run blast for all sequences in the file and return the result in form of elements_file array list
		// so calculate all hits and get top 8
		// CHECK FOR HIT is TM or NOT?
		SimapAccessWebService simap;
		try {
			simap = new SimapAccessWebService();
			
			System.out.print("\nFILE NAME " + parsed_file.get(1).file_name);
			for (int i =0; i<= parsed_file.size()-1; i++){
				if (parsed_file.get(i).tm_no>0){			//for all the tm proteins in file run blast
					queryseqcount++;
					String md5=simap.computeMD5(parsed_file.get(i).seq);
					simap.setMd5(md5);
					System.out.print("\nsequence id " + parsed_file.get(i).id );
					//System.out.print("\nSequence query \n" + parsed_file.get(i).seq + "\n");
					simap.setMax_evalue(10e-25);		//check params
					simap.setMax_number_hits(8);
					simap.setMin_swscore(120);
					
					simap.alignments(true);
					simap.sequences(true);
					simap.excludeTaxon(10239);	//exclude viral hits

					ArrayList<HitSet> result = ResultParser.parseResult(simap.getHitsXML());

					if (!result.isEmpty()){	//to check if any hit exists or not

						int result_iterator = 1;	//1 because the first hit always self ;)

						while(result_iterator<result.size()-1){

							HitSet second=(HitSet) result.get(result_iterator);
							/*

							System.out.println("\n"+"E-Value\n"+second.getHitAlignment().getEvalue());
							System.out.println("\n Sequence Hit: \n"+second.getHitData().getSequence() +"\n");
							System.out.println("\n Sequence Query: \n"+parsed_file.get(i).seq +"\n");
							System.out.println("\n ALLIGNMENT QUERY:\n"+second.getHitAlignment().getAlignment_query() +"\n");
							System.out.println("\n ALLIGNMENT Hit:\n"+second.getHitAlignment().getAlignment_hit() +"\n");
							System.out.println("\n HitStart:\n"+second.getHitAlignment().getHit_start() +"\n");
							System.out.println("\n Hit:\n"+second.getHitAlignment().getHit_stop() +"\n");
							 */
							int qStart = second.getHitAlignment().getQuery_start() -1;
							int qStop = second.getHitAlignment().getQuery_stop() -1;
							
							String alnQ = second.getHitAlignment().getAlignment_query();

							Check_hit_tm obj = new Check_hit_tm();
							boolean pho = obj.run(second.getHitData().getSequence());
							if(pho){	//if hit is a TM then
								
								if (CheckTmCoverage(parsed_file.get(i),qStart,qStop)){	//check if atleast 50% of tm query coverage
									
									hitstatObj = obj.getcounts(second.getHitData().getSequence(),hitstatObj);//was in if above
									hitstatObj.hits_no++;
									
									Elements_file element = getCoverageProfile(parsed_file.get(i),alnQ,qStart,qStop);
									parsed_file.set(i, element);
									//System.out.print("HitId: " + second.getHitData().getMd5()+"\n");
								}
							}
							result_iterator ++;
							// a hit consists out of alignment data and hit data
						}
						//System.out.println("\n*********************************************************************************\n");
					}
				}
			}
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return parsed_file;
	}

	private Elements_file getCoverageProfile(Elements_file e, String aln,int qst, int qstp){
		int j = qst;
		int scorep1 = e.scorep;
		int scoren1 = e.scoren;

		for (int i =0; i<=aln.length()-1 && j<=qstp; i++){

			if (aln.charAt(i)!=e.seq.charAt(j)){
				//for a tm segment---find if a tm segment using j and then while the tm segment ends, check alignment finally it makes sense
				//if j in range of tm segment then while the tm segment ends check each position and score gap and mm:)
				if(inTmOrNot(j,e.tm_no,e.start,e.end)){
					int x = 0;
					int range = e.end[rangeIndex] - e.start[rangeIndex];
					while (x<=range && i<= aln.length()-1 ){// traverse the whole TM to check for changes
						if (aln.charAt(i)!=e.seq.charAt(j)){
							if (e.seq.charAt(j)=='K' || e.seq.charAt(j)=='R' || e.seq.charAt(j)=='H'){
								if (aln.charAt(i) == '-'){
									scorep1 = scorep1 + 2;
								}
								else if(aln.charAt(i) == 'R' || aln.charAt(i) == 'H'||aln.charAt(i) == 'K'){
								//scorep1 = scorep1 + 1;	do nothing because positive residue replaced positive residue
								}
								else{
									scorep1 = scorep1 + 1;
								}
							}

							else if (e.seq.charAt(j)=='D' || e.seq.charAt(j)=='E'){
								if (aln.charAt(i) == '-'){
									scoren1 = scoren1 + 2;
								}
								else if (aln.charAt(i) == 'E' || aln.charAt(i) == 'D'){
									//scoren1 = scoren1 + 1; Do nothing because -ive residue replaced negative
								}
								else
								scoren1 = scoren1 + 1;
							}
						}
						x++;
						j++;
						i++;
					}
				}

			}
			j++;
		}
		e.scorep = scorep1;
		e.scoren = scoren1;
	
		return e;
	}

	private boolean inTmOrNot(int pos,int no, int[] st,int[] end){
		boolean exists = false;
		for (int i =0; i<=no; i++){
			if (pos>= st[i] && pos<=end[i]){
				exists = true;
				rangeIndex = i;
				break;
			}
		}
		return exists;
	}

	private boolean CheckTmCoverage(Elements_file e, int qst, int qstp){	//checks if atleast half of query tm are covered

		int count = 0;
		for (int i =0; i<=e.tm_no-1;i++){
			if (e.start[i] >= qst && e.end[i]<=qstp){
				count ++;
			}
		}
		boolean coverage = false;
		if (count >= e.tm_no/2){
			coverage =true;
		}
		return coverage;
	}
	
	
	private int findGinTM(String aln, int i, int j) { //i tmstart in algn and j is theoretical stop in algn
		// TODO Auto-generated method stub
		String subaln = aln.substring(i, j);
		int count = 0; 
		count= StringUtils.countMatches(subaln, "-");
		if (count>0){
			j = j + count;
			subaln = aln.substring(i, j);
		}
		// now we have original end of tm in algn as j and the segment in subaln
		return 0;
	}

}
