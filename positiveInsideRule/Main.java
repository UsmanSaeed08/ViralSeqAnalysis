package positiveInsideRule;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;

//import org.apache.xpath.operations.String;

import utils.DBAdaptor;

public class Main {

	/**
	 * @param args
	 * @throws SQLException 
	 */

	//static ArrayList <Elements_file> human;

	static HitStats obj_hitstat = new HitStats();
	static int i ;
	private static BufferedWriter wr;
	private static BufferedWriter wrfasta;
	private static long InsideLen = 0;
	private static long TmLen = 0;
	private static long OutsideLen = 0;

	private static final Connection CAMPS_CONNECTION = DBAdaptor.getConnection("CAMPS4");

	private static void Main2(){
		System.out.print("\nWelcome to Main2 \n");
		try{
			final Connection CAMPS_CONNECTION1 = DBAdaptor.getConnection("CAMPS4");
			final Connection CAMPS_CONNECTION2 = DBAdaptor.getConnection("CAMPS4");
			final Connection CAMPS_CONNECTION3 = DBAdaptor.getConnection("CAMPS4");

			HitStats_holder holder = new HitStats_holder();


			Statement stm2 = CAMPS_CONNECTION.createStatement();
			stm2.setFetchSize(Integer.MIN_VALUE);
			ResultSet rs2 = stm2.executeQuery("SELECT sequenceid FROM proteins2 WHERE taxonomyid = 9606"); //9606 humans
			//ResultSet rs2 = stm2.executeQuery("SELECT sequenceid FROM proteins2 WHERE taxonomyid = 10090"); //10090 Musmusculus
			//ResultSet rs2 = stm2.executeQuery("SELECT sequenceid FROM proteins2 WHERE taxonomyid = 6239"); //6239 C.Elegans
			//ResultSet rs2 = stm2.executeQuery("SELECT sequenceid FROM proteins2 WHERE taxonomyid = 3702"); //3702 Arabidiopsis
			//ResultSet rs2 = stm2.executeQuery("SELECT sequenceid FROM proteins2 WHERE taxonomyid = 7227"); //7227 Drosophila
			System.out.print("\n Done with Connection stuff \n");
			int numberOfprots = 0;
			while(rs2.next()) {
				numberOfprots++;
				int sequenceid = rs2.getInt("sequenceid");
				//PreparedStatement pstm = CAMPS_CONNECTION.prepareStatement("SELECT begin,end FROM elements WHERE sequenceid=?");
				PreparedStatement pstm = CAMPS_CONNECTION1.prepareStatement("SELECT begin,end FROM tms WHERE sequenceid=?");
				PreparedStatement pstm_seq = CAMPS_CONNECTION2.prepareStatement("SELECT sequence FROM sequences2 WHERE sequenceid=?");
				PreparedStatement pstm_topology = CAMPS_CONNECTION3.prepareStatement("SELECT n_term,c_term FROM topology WHERE sequenceid=?");

				pstm.setInt(1, sequenceid);
				pstm_seq.setInt(1,sequenceid);
				pstm_topology.setInt(1,sequenceid);


				ResultSet rs_pstm_seq = pstm_seq.executeQuery();

				String seq ="";
				
				// GET SEQUENCE
				while(rs_pstm_seq.next()){
					seq = rs_pstm_seq.getString("sequence");				
				}
				rs_pstm_seq.close();
				pstm_seq.close();

				ResultSet rs_pstm = pstm.executeQuery();
				Elements_file e = new Elements_file();

				// GET TMs BEGIN AND END
				while(rs_pstm.next()){
					//Elements_file e = new Elements_file();
					//System.out.print(rs_pstm.getString("md5")+"\n");
					e.id = String.valueOf(sequenceid);
					e.start[e.tm_no] = rs_pstm.getInt("begin");
					e.end[e.tm_no] = rs_pstm.getInt("end");
					e.seq = seq;
					//String segment = seq.substring(e.start[e.tm_no-1], e.end[e.tm_no-1]);	//-1 because starts from zero
					//e.tm_segments.add(segment);
					e.tm_no = e.tm_no + 1;
					//obj_hitstat = obj.getcounts(seq, obj_hitstat, e.start, e.end, e.tm_no);
				}
				rs_pstm.close();
				pstm.close();

				ResultSet rs_pstm_topology = pstm_topology.executeQuery();

				// GET TOPOLOGY
				while(rs_pstm_topology.next()){
					if (rs_pstm_topology.getString(1).contains("in")){
						e.nterm = 1;
					}
					else if (rs_pstm_topology.getString(1).contains("out")){
						e.nterm = 0;
					}

					if (rs_pstm_topology.getString(2).contains("in")){
						e.cterm = 1;
					}
					else if (rs_pstm_topology.getString(2).contains("out")){
						e.cterm = 0;
					}
				}
				rs_pstm_topology.close();
				pstm_topology.close();

				// get all the segments
				ArrayList<String> Segments = new ArrayList<String>();
				// j is the current tm_no
				String tempSeq = e.seq;
				for(int j =0; j<=e.tm_no-1;j++){
					
					//tmSegments.add(e.seq.substring(e.start[j], e.end[j]));
					int st = e.start[j]; 
					int end = e.end[j];
					//if (st - end <5)
						//continue;
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
					// edit the sequence for the parts which have been attained

					//tempSeq = tempSeq.substring(end-1);
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
						continue;
					}
					else if(e.cterm == 0 && lastEntry != 3){
						// bad
						// no use of the whole protein here.. so it is useless to count all the hit stats, so just continue
						continue;
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
						continue;
					}
					else if(e.cterm == 0 && lastEntry != 3){
						// bad
						// no use of the whole protein here.. so it is useless to count all the hit stats, so just continue
						continue;
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
					//thisprotein.obj_hitstats_outsdide = Check_hit_tm.count(TmSegments.get(x), thisprotein.obj_hitstats_outsdide);
					//thisprotein.obj_hitstats_outsdide.NumberOfthisSegment = thisprotein.obj_hitstats_outsdide.NumberOfthisSegment + 1 ;
					
				}
				for (int x = 0;x <= OutsideSegments.size()-1; x++){
					thisprotein.obj_hitstats_outsdide = Check_hit_tm.count(OutsideSegments.get(x), thisprotein.obj_hitstats_outsdide);
					thisprotein.obj_hitstats_outsdide.NumberOfthisSegment = thisprotein.obj_hitstats_outsdide.NumberOfthisSegment + 1 ;
				}
				
				/*
				// get percentage of residues 
				// i.e divide frequency of residue by total residues
				thisprotein.obj_hitstats_inside = Check_hit_tm.makePercentage(thisprotein.obj_hitstats_inside);
				thisprotein.obj_hitstats = Check_hit_tm.makePercentage(thisprotein.obj_hitstats);
				thisprotein.obj_hitstats_outsdide = Check_hit_tm.makePercentage(thisprotein.obj_hitstats_outsdide);
				// now taking percentage.. testing to see what happens if we normalize by the number of such segments
				thisprotein.obj_hitstats_inside = Check_hit_tm.makePercentage2(thisprotein.obj_hitstats_inside,InsideSegments.size());
				thisprotein.obj_hitstats = Check_hit_tm.makePercentage2(thisprotein.obj_hitstats,TmSegments.size());
				// since i combine tm with outside above, so its should be divided by more
				thisprotein.obj_hitstats_outsdide = Check_hit_tm.makePercentage2(thisprotein.obj_hitstats_outsdide,OutsideSegments.size());
				*/
				
				// now the array has percentages and indivisual integers have count
				// now add thisprotein to holder because otherwise all the ercentages and counts would
				// get mixed up
				holder.obj_hitstats_inside = Check_hit_tm.SumUp(holder.obj_hitstats_inside,thisprotein.obj_hitstats_inside);
				holder.obj_hitstats = Check_hit_tm.SumUp(holder.obj_hitstats,thisprotein.obj_hitstats);
				holder.obj_hitstats_outsdide = Check_hit_tm.SumUp(holder.obj_hitstats_outsdide,thisprotein.obj_hitstats_outsdide);

				System.out.print(".");
				if (i % 100 == 0){
					System.out.print("\n" + i + "\n");
				}
				i++;
				//write_toFile(e);
				//write_FastatoFile(e);
			}
			//numberOfprots

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

			System.out.print("No of Prots "+numberOfprots+"\n");

			stm2.close();
			rs2.close();

			CAMPS_CONNECTION1.close();
			CAMPS_CONNECTION2.close();
			CAMPS_CONNECTION3.close();
			CAMPS_CONNECTION.close();


		}
		catch(Exception e){
			e.printStackTrace();
		}

	}

	public static void main(String[] args) throws SQLException, IOException {
		// TODO Auto-generated method stub

		Main2(); // made to test the positive inside rule all over again

		//String path = "F:/virus_project/phobius_results.individualproteins";
		//String path = "F:/virus_project/1/tm";
		//String c = "check";
		//Fetch_Directory obj = new Fetch_Directory(path);

		/*
		1. get sp from virus tm..
		2. get the sp segments
		3. aa composition --> sys.out
		 */

		// ******************************* READ TODB XML FILE START ***************************
		/*
		ReadXMLFile_TopDb obj = new ReadXMLFile_TopDb();
		HitStats object = new HitStats();
		object = obj.run(object);
		System.out.print("\n Stats are \n \n");
		for (int i =0; i<=19; i++){
			System.out.print(object.aa[i]+"\t");
		}
		//for (int i =0; i<=ReadXMLFile_TopDb.virus-1; i++){
			//System.out.print("\n"+ReadXMLFile_TopDb.names.get(i)+"\n");
		//}
		System.out.print("\n Viral Proteins Analyzed : "+ReadXMLFile_TopDb.virus_TMProteins+"\n");
		System.out.print("\n Viral TM Regions Analyzed : "+ReadXMLFile_TopDb.virus_TMRegions+"\n");

		System.out.print("\n Cellular Proteins Analyzed : "+ReadXMLFile_TopDb.cellular_Proteins+"\n");
		System.out.print("\n Cellular TM Regions Analyzed : "+ReadXMLFile_TopDb.cellular_TMRegions+"\n");

		System.out.print("\n Human: " + ReadXMLFile_TopDb.h);
		System.out.print("\n Mouse: " + ReadXMLFile_TopDb.m);
		System.out.print("\n EColi: " + ReadXMLFile_TopDb.E);
		System.out.print("\n Drosophila: " + ReadXMLFile_TopDb.D);
		System.out.print("\n Yeast: " + ReadXMLFile_TopDb.S);
		System.out.print("\n CElegans: " + ReadXMLFile_TopDb.C);
		System.out.print("\n Arabidiopsis: " + ReadXMLFile_TopDb.A);

		 */
		//System.out.print("\n Cellular Proteins Analyzed : "+ReadXMLFile_TopDb.cellular+"\n");

		// ******************************* READ TODB XML FILE END ***************************


		//************************************************ EXTRACT HUMAN AND OTHER DATA START *********************************

		/*

		final Connection CAMPS_CONNECTION = DBAdaptor.getConnection("CAMPS4");

		//wr = new BufferedWriter(new FileWriter("F:/drosophila_tm.txt"));
		//wrfasta = new BufferedWriter(new FileWriter("F:/helix interaction graph project/TMproteins/DrosophilaTm.txt"));

		i = 0;
		Statement stm2 = CAMPS_CONNECTION.createStatement();
		//Statement stm3 = CAMPS_THREE_CONNECTION.createStatement();
		//stm2.setFetchSize(Integer.MIN_VALUE);
		//ResultSet rs2 = stm2.executeQuery("SELECT sequenceid FROM proteins WHERE taxonomyid = 9606"); //9606 humans
		//ResultSet rs2 = stm2.executeQuery("SELECT sequenceid FROM proteins2 WHERE taxonomyid = 10090"); //10090 Musmusculus
		//ResultSet rs2 = stm2.executeQuery("SELECT sequenceid FROM proteins2 WHERE taxonomyid = 6239"); //6239 C.Elegans
		//ResultSet rs2 = stm2.executeQuery("SELECT sequenceid FROM proteins2 WHERE taxonomyid = 3702"); //3702 Arabidiopsis
		ResultSet rs2 = stm2.executeQuery("SELECT sequenceid FROM proteins2 WHERE taxonomyid = 7227"); //7227 Drosophila


		//ResultSet rs3 = stm3.executeQuery("SELECT sequences_sequenceid FROM camps2uniprot WHERE taxonomyid = 4932"); //4932 Yeast
		//ResultSet rs3 = stm3.executeQuery("SELECT sequences_sequenceid FROM camps2uniprot WHERE taxonomyid = 562"); //562 Yeast
		//human = new Elements_file[size];

		//human = new ArrayList <Elements_file> ();


		System.out.print("\n Done with Connection stuff \n");
		while(rs2.next()) {
			int sequenceid = rs2.getInt("sequenceid");


			//PreparedStatement pstm = CAMPS_CONNECTION.prepareStatement("SELECT begin,end FROM elements WHERE sequenceid=?");
			PreparedStatement pstm = CAMPS_CONNECTION.prepareStatement("SELECT begin,end FROM tms WHERE sequenceid=?");
			PreparedStatement pstm_seq = CAMPS_CONNECTION.prepareStatement("SELECT sequence FROM sequences2 WHERE sequenceid=?");

			pstm.setInt(1, sequenceid);
			pstm_seq.setInt(1,sequenceid);

			ResultSet rs_pstm = pstm.executeQuery();
			ResultSet rs_pstm_seq = pstm_seq.executeQuery();

			String seq ="";
			Check_hit_tm obj = new Check_hit_tm();

			while(rs_pstm_seq.next()){
				seq = rs_pstm_seq.getString("sequence");				
			}
			/*
		 * changes here to just get the sequences for the run of hblits and predict contacts





			//if(seq != ""){

			//IF YOU WANT TO GET THE SP AND SEGMENTS OF TM INFO THEN DO BELOW

			Elements_file e = new Elements_file();
			while(rs_pstm.next()){
				//Elements_file e = new Elements_file();
				//System.out.print(rs_pstm.getString("md5")+"\n");
				e.id = String.valueOf(sequenceid);

				e.start[e.tm_no] = rs_pstm.getInt("begin");
				e.end[e.tm_no] = rs_pstm.getInt("end");
				e.seq = seq;
				//String segment = seq.substring(e.start[e.tm_no-1], e.end[e.tm_no-1]);	//-1 because starts from zero
				//e.tm_segments.add(segment);

				e.tm_no = e.tm_no + 1;

				//obj_hitstat = obj.getcounts(seq, obj_hitstat, e.start, e.end, e.tm_no);

				System.out.print(".");
				if (i % 100 == 0){
					System.out.print("\n" + i + "\n");
				}
				//now get the sp segment using begin and end
				i++;
			}
			//write_toFile(e);
			write_FastatoFile(e);

		}
		//wr.close();
		//wrfasta.close();

		//}

		 */
		//************************************************ EXTRACT HUMAN AND OTHER DATA START *********************************


		//String path = "D:/virus_project/phobius_results.individualproteins";
		//String path = "D:/virus_project/1/tm";
		//Fetch_Directory obj = new Fetch_Directory(path);

		// so now have to get for all other organisms, starting with human we dont see the viral inputs!
		//steps.... 1. get seq from camps.... 2. use sequences to get phobius from simap.... 3. calc stats

		/*
		Check_hit_tm obj = new Check_hit_tm();
		String seq= "MEKWYLMTVVVLIGLTVRWTVSLNSYSGAGKPPMFGDYEAQRHWQEITFNLPVKQWYFNSSDNNLQYWGLDYPPLTAYHSLLCAYVAKFINPDWIALHTSRGYESQAHKLFMRTTVLIADLLIYIPAVVLYCCCLKEISTKKKIANALCILLYPGLILIDYGHFQNIYNSVSLGFALWGVLGISCDCDLLGSLAFCLAINYKQMELYHALPffcfllgkcfkkglkgkgfvllvklACIVVASFVLCWLPFFTEREQTLQVLRRLFPVDRGLFEDKVANIWCSFNVFLKIKDILPRHIQLIMSFCSTFLSLLPACIKLILQPSSKGFKFTLVSCALSFFLFSFQVHEKsillvslpvclvlsEIPFMSTWFLLVSTFSMlplllkdellMPSVVTTMAFFIACVTSFSIFEKTSEEELQLKSFSISVRKYLPCFTFLSRIIQYLFLISVItmvlltlmtvtlDPPQKLPDLFSVLVCFVSCLNFLFFLVYFNIIIMWDSKSGRNQKKIS";
		obj.run(seq);
		System.out.print("Length " + seq.length() + "\n");
		for(int i=0; i<=obj.tmsCount-1;i++){
			System.out.print("st " + obj.start_[i] + "\t");
			System.out.print("end " + obj.end_[i] + "\n");
		}

		final Connection CAMPS_CONNECTION = DBAdaptor.getConnection("CAMPS4");
		final Connection CAMPS_THREE_CONNECTION = DBAdaptor.getConnection("CAMPS3");

		i = 0;
		Statement stm2 = CAMPS_CONNECTION.createStatement();
		//Statement stm3 = CAMPS_THREE_CONNECTION.createStatement();
		//stm2.setFetchSize(Integer.MIN_VALUE);
		ResultSet rs2 = stm2.executeQuery("SELECT sequenceid FROM proteins WHERE taxonomyid = 9606"); //9606 humans
		//ResultSet rs2 = stm2.executeQuery("SELECT sequenceid FROM proteins WHERE taxonomyid = 10090"); //10090 Musmusculus
		//ResultSet rs2 = stm2.executeQuery("SELECT sequenceid FROM proteins WHERE taxonomyid = 6239"); //6239 C.Elegans
		//ResultSet rs2 = stm2.executeQuery("SELECT sequenceid FROM proteins WHERE taxonomyid = 3702"); //3702 Arabidiopsis
		//ResultSet rs2 = stm2.executeQuery("SELECT sequenceid FROM proteins WHERE taxonomyid = 7227"); //7227 Drosophila


		//ResultSet rs3 = stm3.executeQuery("SELECT sequences_sequenceid FROM camps2uniprot WHERE taxonomyid = 4932"); //4932 Yeast
		//ResultSet rs3 = stm3.executeQuery("SELECT sequences_sequenceid FROM camps2uniprot WHERE taxonomyid = 562"); //562 Yeast
		//human = new Elements_file[size];

		//human = new ArrayList <Elements_file> ();


		System.out.print("\n Done with Connection stuff \n");
		while(rs2.next()) {
			int sequenceid = rs2.getInt("sequenceid");


			PreparedStatement pstm = CAMPS_CONNECTION.prepareStatement("SELECT md5,sequence FROM sequences WHERE sequenceid=?");
			pstm.setInt(1, sequenceid); 
			ResultSet rs_pstm = pstm.executeQuery();
			while(rs_pstm.next()){
				Elements_file e = new Elements_file();
				//System.out.print(rs_pstm.getString("md5")+"\n");
				e.id = String.valueOf(sequenceid);
				e.md5 = rs_pstm.getString("md5");
				e.seq = rs_pstm.getString("sequence");
				i++;
				get_simap_tm(e.seq);

			}
		}
		//HitStats h = get_simap_tm(size);
		//size = human.size();
		System.out.print("\n Extraction of tm sections done \n");
		//HitStats obj_hitstat = new HitStats();
		//get_stat(size,obj_hitstat);
		//System.out.print("\n Stats gathered \n");
		// Now to make the charts...
		 */

		// ************************** TO PRINT THE STATS *****************************************************
		/*
		System.out.print("\n TOTAL SP no: " + i + "\n");
		for(int i =0; i<= obj_hitstat.aa_seq.length()-1;i++){
			System.out.print("\n"+ obj_hitstat.aa_seq.charAt(i) + " --- " + obj_hitstat.aa[i]);
		}
		System.out.print("\n"+"no of SP: " + obj_hitstat.hits_no);

		//Fetch_Directory obj = new Fetch_Directory();
		//obj.makechart(obj_hitstat);
		 */
	}


	// to write all details to file
	private static void write_toFile(Elements_file e) throws IOException {
		// TODO Auto-generated method stub
		wr.write("ID	" + e.id);
		wr.newLine();
		wr.write("SEQ	" + e.seq);
		wr.newLine();
		for (int i = 0; i<= e.tm_no-1; i++){
			wr.write("FT	TRANSMEME	" + e.start[i]+ "	"+ e.end[i]);
			wr.newLine();
		}
		wr.write("//");
		wr.newLine();
	}
	//to write only fasta sequences and ids to file
	private static void write_FastatoFile(Elements_file e) throws IOException {
		// TODO Auto-generated method stub
		try{
			BufferedWriter writer;
			writer = new BufferedWriter(new FileWriter("F:/helix interaction graph project/TMproteins/DrosophilaSeq/"+e.id+".fasta"));
			writer.write(">" + e.id);
			writer.newLine();
			writer.write(e.seq);
			writer.newLine();
			writer.close();
		}
		catch(Exception ex){
			ex.printStackTrace();
		}

	}

	private static void get_simap_tm(String sq) {
		// TODO Auto-generated method stub
		Check_hit_tm obj = new Check_hit_tm();

		if(obj.run(sq)){
			Elements_file e = new Elements_file();
			e.start = obj.start_;
			e.end = obj.end_;
			e.tm_no = obj.tmsCount;

			obj_hitstat = obj.getcounts(sq, obj_hitstat, e.start, e.end, e.tm_no);

			System.out.print(".");
			if (i % 100 == 0){
				System.out.print("\n" + i + "\n");
			}
		}

	}

}
