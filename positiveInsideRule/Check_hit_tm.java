package positiveInsideRule;



import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.logging.Logger;
import java.util.zip.GZIPInputStream;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;

import mips.gsf.de.simapclient.client.SimapAccessWebService;
import mips.gsf.de.simapfeatureclient.client.FeatureClient;

import org.w3c.dom.Document;
import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

public class Check_hit_tm {
	
	public boolean phobius = false;

	public int start_[] = new int[200];
	public int end_[] = new int[200];
	int tmsCount = 0;

	public boolean run(String seq) { // seq is the hit sequence for which it has to be checked if the tm exists?

		//<tns:details/></feature><feature description="Signal peptide" evidence="PHOBIUS" modelid="PHOBIUS_SP" modelref="PHOBIUS" 
		
		int searchedCount = 0;
		int tm1count = 0;
		

			searchedCount++;
			boolean p = false;

			FeatureClient f = null;
			SimapAccessWebService simap = null;
			Document doc = null;

			try {
				f = new FeatureClient("http://ws.csb.univie.ac.at/simapwebservice/services/SimapService");
				simap = new SimapAccessWebService();
				String md5 = simap.computeMD5(seq);
				DocumentBuilderFactory dbFactory = DocumentBuilderFactory.newInstance();
				DocumentBuilder dBuilder = dbFactory.newDocumentBuilder();
				doc = dBuilder.parse(new ByteArrayInputStream(f.getFeaturesXML(md5).getBytes()));



			NodeList featureNodes = doc.getElementsByTagName("feature");
			for (int i = 0; i < featureNodes.getLength(); i++) {
				Node featureNode = featureNodes.item(i);
				NamedNodeMap featureNodeAttributes = featureNode.getAttributes();
				phobius = false;
				for (int j = 0; j < featureNodeAttributes.getLength(); j++) {	//check if phobius prediction
					Node featureNodeAttribute = featureNodeAttributes.item(j);
					
					//if (featureNodeAttribute.getNodeName().equals("modelid") && featureNodeAttribute.getNodeValue().equals("PHOBIUS_TM")) {
					if (featureNodeAttribute.getNodeName().equals("modelid") && featureNodeAttribute.getNodeValue().equals("PHOBIUS_SP")) {
						phobius = true;
						p = true;
						//System.out.print(f.getFeaturesXML(md5));
						break;
						
					}
				}
				/*
				if (phobius){
					System.out.print(f.getFeaturesXML(md5));
				}*/

				if (phobius) {	//only if phobius prediction exists for the hit
					NodeList locationNodes = featureNode.getChildNodes();
					
					int start = -1;
					
					int z =-1;
					int end = -1;
					for (int j = 0; j < locationNodes.getLength(); j++) {
						Node locationNode = locationNodes.item(j);
						if (locationNode.getNodeName().equals("location")) {
							tmsCount++;
							z++;
							NodeList posNodes = locationNode.getChildNodes();
							for (int k = 0; k < posNodes.getLength(); k++) {
								Node posNode = posNodes.item(k);
								if (posNode.getNodeName().equals("begin")) {
									NamedNodeMap beginNodeAttributes = posNode.getAttributes();
									for (int l = 0; l < beginNodeAttributes.getLength(); l++) {
										Node beginNodeAttribute = beginNodeAttributes.item(l);
										if (beginNodeAttribute.getNodeName().equals("position")) {
//											System.out.println("begin: " + beginNodeAttribute.getNodeValue());
											start = Integer.parseInt(beginNodeAttribute.getNodeValue()) - 1;
											start_[z] = start;
										}
									}
								}
								if (posNode.getNodeName().equals("end")) {
									NamedNodeMap endNodeAttributes = posNode.getAttributes();
									for (int l = 0; l < endNodeAttributes.getLength(); l++) {
										Node endNodeAttribute = endNodeAttributes.item(l);
										if (endNodeAttribute.getNodeName().equals("position")) {
//											System.out.println("end: " + endNodeAttribute.getNodeValue());
											end = Integer.parseInt(endNodeAttribute.getNodeValue()) - 1;
											end_[z] = end;
										}
									}
								}
							}
						}
					}
					int length = end - start + 1;
					//if (tmsCount == 1) {
					//	tm1count++;
					//}
				}
			//	if (phobius){	//again...if true and till now the start and end have been calculated
								// then we can not get a sequence profile
					
					
				//}
			//	return phobius;
				
			}
			
			
			//System.out.println("\n done with finding hit tm \n");	//by here phobius gets false!!!
			//return phobius;
			} catch (Exception ex) {
				ex.printStackTrace();

			}
		
			return p;
 
		//log.info("1 TM proteins: " + tm1count);

	}
	public HitStats getcounts(String seq, HitStats oj){
			// if the hit has predictions then get stats for hit tm segments
			for (int i =0; i<= tmsCount-1; i++){
				String tm = seq.substring(start_[i], end_[i]);	//geting sequence
				oj = count(tm, oj);
				//element.tm_segments.add(tm);
			}
		return oj;
	}
	public HitStats getcounts(String seq, HitStats oj, int st[], int end[], int tmNo){
		// if the hit has predictions then get stats for hit tm segments
		for (int i =0; i<= tmNo-1; i++){
			
			if(st[i] < seq.length()-1 && end[i] < seq.length()-1){
			String tm = seq.substring(st[i]-1, end[i]-1);	//geting sequence
			oj = count(tm, oj);
			}
			else 
				System.out.print("\n" + seq +"\n"+ "St " + st[i] + "------- En " + end[i] +"\n");
			//element.tm_segments.add(tm);
		}
	return oj;
}
	public HitStats_holder getcounts_rule_one(String seq, HitStats_holder oj, int st[], int end[], int tmNo, Elements_file p){
		// if the hit has predictions then get stats for hit tm segments
//		if (tmNo == 1){
			if(seq.length()-1 > st[0]-1){
				
			
			String inside = seq.substring(0, st[0]-1);
			String tm = seq.substring(st[0]-1, end[0]-1);
			String outside = seq.substring(end[0]-1, seq.length()-1);
			
			oj.obj_hitstats_inside = count(inside, oj.obj_hitstats_inside);
			oj.obj_hitstats_outsdide = count(outside, oj.obj_hitstats_outsdide);
			oj.obj_hitstats = count(tm, oj.obj_hitstats);
			}
			else
				System.out.print("\n \n REPORT st exception\n \n " + seq + "\n" + st[0] + "\n" + p.file_name + "\n");
	//	}
		
	return oj;
}
	public HitStats_holder getcounts_rule(String seq, HitStats_holder oj, int st[], int end[], int tmNo){
		// if the hit has predictions then get stats for hit tm segments
		boolean flag = true;
		String inside = "";
		String tm = "";
		String outside = "";
		
		for (int i =0; i<= tmNo-1; i++){
			
			if (flag){
				if (i==0){
					inside = seq.substring(0, st[0]-1);
					tm = seq.substring(st[0]-1, end[0]-1);
					
					oj.obj_hitstats_inside = count(inside, oj.obj_hitstats_inside);
					oj.obj_hitstats = count(tm, oj.obj_hitstats);
					
				}
				else{
					inside = seq.substring(end[i-1], st[i]-1);
					tm = seq.substring(st[i]-1, end[i]-1);
					
					oj.obj_hitstats_inside = count(inside, oj.obj_hitstats_inside);
					oj.obj_hitstats = count(tm, oj.obj_hitstats);
				}
				flag = false;
				
			}
			else if(!flag){
				
				tm = seq.substring(st[i], end[i]);
				outside = seq.substring(end[i-1]+1, st[i]);
				
				oj.obj_hitstats_outsdide = count(outside, oj.obj_hitstats_outsdide);
				oj.obj_hitstats = count(tm, oj.obj_hitstats);
				flag = true;
			}
			if (i == tmNo-1){
				inside = seq.substring(end[i]+1, seq.length()-1);
				oj.obj_hitstats_inside = count(inside, oj.obj_hitstats_inside);
			}
			
		}
	return oj;
}
	public static HitStats makePercentage(HitStats obj){
		for(int i = 0; i<=obj.aa_seq.length()-1;i++){
			if(obj.totalResidues>0){
			obj.aa[i] = obj.aa[i]/obj.totalResidues;
			//obj.aa[i] = obj.aa[i] * 100;
			obj.aa[i] = obj.aa[i];
			
			obj.hydrophobic = obj.hydrophobic/obj.totalResidues;
			obj.hydrophobic = obj.hydrophobic * 100;
			
			obj.hydrophilic = obj.hydrophilic/obj.totalResidues;
			obj.hydrophilic = obj.hydrophilic * 100;
			
			//
			
			obj.KnR = obj.KnR/obj.totalResidues;
			//obj.KnR = obj.KnR * 100;
			
			}
		}
		return obj;
	}
	public static HitStats makePercentage2(HitStats obj,int segments){
		for(int i = 0; i<=obj.aa_seq.length()-1;i++){
			if(segments>0){
			obj.aa[i] = obj.aa[i]/segments;
			obj.aa[i] = obj.aa[i] * 100;
			
			obj.hydrophobic = obj.hydrophobic/segments;
			obj.hydrophobic = obj.hydrophobic * 100;
			
			obj.hydrophilic = obj.hydrophilic/segments;
			obj.hydrophilic = obj.hydrophilic * 100;
			
			obj.KnR = obj.KnR/segments;
			obj.KnR = obj.KnR * 100;
			}
			else {
				obj.aa[i] = 0; 
			}
		}
		return obj;
	}
	public static HitStats SumUp(HitStats obj,HitStats obj_thisprotein){
		
		obj.totallength = obj.totallength + obj_thisprotein.totallength;
		obj.hydrophobic = obj.hydrophobic + obj_thisprotein.hydrophobic;
		obj.hydrophilic = obj.hydrophilic + obj_thisprotein.hydrophilic;
		obj.KnR = obj.KnR + obj_thisprotein.KnR;
		
		for(int i = 0; i<=obj.aa_seq.length()-1;i++){
			obj.aa[i] = obj.aa[i]+obj_thisprotein.aa[i];
		}
		return obj;
	}
	public static HitStats count(String transmem, HitStats obj){
		obj.totallength = obj.totallength + transmem.length();
		for (int i = 0; i<=transmem.length()-1;i++){
			obj.totalResidues = obj.totalResidues + 1; 
			if (transmem.charAt(i) == 'G'){
				obj.G++;
				obj.aa[0] = obj.aa[0]+1;
				obj.hydrophobic++;
			}
			else if(transmem.charAt(i) == 'A'){
				obj.A++;
				obj.aa[1] = obj.aa[1]+1;
				obj.hydrophobic++;
			}
			else if(transmem.charAt(i) == 'V'){
				obj.V++;
				obj.aa[2] = obj.aa[2]+1;
				obj.hydrophobic++;
			}
			else if(transmem.charAt(i) == 'L'){
				obj.L++;
				obj.aa[3] = obj.aa[3]+1;
				obj.hydrophobic++;
			}
			else if(transmem.charAt(i) == 'I'){
				obj.I++;
				obj.aa[4] = obj.aa[4]+1;
				obj.hydrophobic++;
			}
			else if(transmem.charAt(i) == 'S'){
				obj.S++;
				obj.aa[5] = obj.aa[5]+1;
				obj.hydrophilic++;
			}
			else if(transmem.charAt(i) == 'T'){
				obj.T++;
				obj.aa[6] = obj.aa[6]+1;
				obj.hydrophilic++;
			}
			else if(transmem.charAt(i) == 'C'){
				obj.C++;
				obj.aa[7] = obj.aa[7]+1;
			}
			else if(transmem.charAt(i) == 'M'){
				obj.M++;
				obj.aa[8] = obj.aa[8]+1;
				obj.hydrophobic++;
			}
			else if(transmem.charAt(i) == 'P'){
				obj.P++;
				obj.aa[9] = obj.aa[9]+1;
				obj.hydrophobic++;
			}
			else if(transmem.charAt(i) == 'H'){
				obj.H++;
				obj.aa[10] = obj.aa[10]+1;
				obj.hydrophilic++;
			}
			else if(transmem.charAt(i) == 'R'){
				obj.R++;
				obj.aa[11] = obj.aa[11]+1;
				obj.KnR = obj.KnR + 1;
			}
			else if(transmem.charAt(i) == 'N'){
				obj.N++;
				obj.aa[12] = obj.aa[12]+1;
				obj.hydrophilic++;
			}
			else if(transmem.charAt(i) == 'Q'){
				obj.Q++;
				obj.aa[13] = obj.aa[13]+1;
				obj.hydrophilic++;
			}
			else if(transmem.charAt(i) == 'E'){
				obj.E++;
				obj.aa[14] = obj.aa[14]+1;
			}
			else if(transmem.charAt(i) == 'D'){
				obj.D++;
				obj.aa[15] = obj.aa[15]+1;
			}
			else if(transmem.charAt(i) == 'F'){
				obj.F++;
				obj.aa[16] = obj.aa[16]+1;
				obj.hydrophobic++;
			}
			else if(transmem.charAt(i) == 'W'){
				obj.W++;
				obj.aa[17] = obj.aa[17]+1;
				obj.hydrophobic++;
			}
			else if(transmem.charAt(i) == 'Y'){
				obj.Y++;
				obj.aa[18] = obj.aa[18]+1;
				obj.hydrophilic++;
			}
			else if(transmem.charAt(i) == 'K'){
				obj.K++;
				obj.aa[19] = obj.aa[19]+1;
				obj.KnR = obj.KnR + 1;
			}
		}
		return obj;
	}
	public HitStats getcounts_sp(String seq, HitStats oj,int sp_begin, int sp_end) {
		// TODO Auto-generated method stub
		
			if(sp_begin < seq.length()-1 && sp_end < seq.length()-1){
			String tm = seq.substring(sp_begin-1, sp_end-1);	//geting sequence
			oj = count(tm, oj);
			}
			else 
				System.out.print("\n" + seq +"\n"+ "St " + sp_begin + "------- En " + sp_end +"\n");
			//element.tm_segments.add(tm);
		
	return oj;
		
	}

}
