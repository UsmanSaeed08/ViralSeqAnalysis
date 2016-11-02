package positiveInsideRule;

import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.DocumentBuilder;
import org.w3c.dom.Document;
import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.NodeList;
import org.w3c.dom.Node;
import org.w3c.dom.Element;
import java.io.File;
import java.util.ArrayList;
import java.util.List;


public class ReadXMLFile_TopDb {
	public static int virus =0;
	public static int virus_TMRegions =0;
	public static int virus_TMProteins =0;
	
	public static int cellular =0;
	public static int cellular_TMRegions =0;
	public static int cellular_Proteins =0;
	
	public static int h = 0;
	public static int m = 0;
	public static int A = 0;
	public static int S = 0;
	public static int D = 0;
	public static int E = 0;
	public static int C = 0;
	
	public static List <String> names = new ArrayList<String>();

	public HitStats run(HitStats oj) {
		try {

			File fXmlFile = new File("F:/virus_project/ViralvsCelTM_swissprot/topdb_all.xml");
			DocumentBuilderFactory dbFactory = DocumentBuilderFactory.newInstance();
			DocumentBuilder dBuilder = dbFactory.newDocumentBuilder();
			Document doc = dBuilder.parse(fXmlFile);

			doc.getDocumentElement().normalize();

			System.out.println("Root element :" + doc.getDocumentElement().getNodeName());

			NodeList nList = doc.getElementsByTagName("TOPDB");

			System.out.println("----------------------------");

			for (int temp = 0; temp < nList.getLength(); temp++) {

				Node nNode = nList.item(temp);

				System.out.println("\nCurrent Element :" + nNode.getNodeName());

				if (nNode.getNodeType() == Node.ELEMENT_NODE) {

					Element eElement = (Element) nNode;

					System.out.println("ID : " + eElement.getAttribute("ID"));
					String org = eElement.getElementsByTagName("Organism").item(0).getTextContent();
					System.out.println("Org : " + org);
					org = org.trim();
					org = org.replaceAll(" ", "");
					org = org.replaceAll("\n", "");

					if (org.equals("Homosapiens")||org.equals("Musmusculus")||org.equals("Arabidopsisthaliana")||org.equals("Caenorhabditiselegans")||org.equals("Drosophilamelanogaster")||org.equals("Escherichiacoli")||org.equals("Saccharomycescerevisiae")){
					//if (org.contains("virus")||org.contains("Virus")){
						
						if (org.equals("Homosapiens")){
							h++; 
						}
						else if(org.equals("Musmusculus")){
							m++;
						}
						else if(org.equals("Arabidopsisthaliana")){
							A++;
						}
						else if(org.equals("Saccharomycescerevisiae")){
							S++;				
						}
						else if(org.equals("Escherichiacoli")){
							E++;
						}
						else if(org.equals("Drosophilamelanogaster")){
							D++;
						}
						else if(org.equals("Caenorhabditiselegans")){
							C++;
						}
						
					ReadXMLFile_TopDb.virus++;
					ReadXMLFile_TopDb.names.add(org);
					
					//ReadXMLFile_TopDb.cellular++;
						
					String seq = eElement.getElementsByTagName("Sequence").item(0).getTextContent();
					seq = seq.trim();
					seq =seq.replaceAll(" ", "");
					seq =seq.replaceAll("\n", "");

					System.out.println("Seq : " + seq);

					NodeList n =nNode.getChildNodes();
					for (int i=0;i<=n.getLength()-1;i++){
						Node mynode = n.item(i);
						if (mynode.getNodeName().equals("Topology")){
							ReadXMLFile_TopDb.virus_TMProteins++;
							ReadXMLFile_TopDb.cellular_Proteins++;
							
							NodeList n2 = mynode.getChildNodes();
							for (int j =0; j<=n2.getLength()-1;j++){
								Node mynode2 = n2.item(j);
								if (mynode2.getNodeName().equals("Regions")){
									NodeList n3 = mynode2.getChildNodes();
									for (int z =0; z<=n3.getLength()-1;z++){
										Node reg = n3.item(z);
										if (reg.getNodeType() == Node.ELEMENT_NODE) {
											Element ReElement = (Element) reg;
											//System.out.println("Begin : " + ReElement.getAttribute("Begin") + "\t End : " + ReElement.getAttribute("End") + "\t Loc : " + ReElement.getAttribute("Loc"));
											String begin = ReElement.getAttribute("Begin");
											String end = ReElement.getAttribute("End");
											String loc = ReElement.getAttribute("Loc");
											
											begin = begin.trim();
											end = end.trim();
											loc = loc.trim();
											
											if (loc.equals("Membrane")){	// just to avoid in consistency
												ReadXMLFile_TopDb.virus_TMRegions++;
												ReadXMLFile_TopDb.cellular_TMRegions++;
												
												Integer bg = Integer.parseInt(begin);
												Integer ed = Integer.parseInt(end);
												String tm = seq.substring(bg, ed);
												Check_hit_tm hitTMCheck = new Check_hit_tm();
												oj = hitTMCheck.count(tm, oj);
											}
											
										}
									}
								}


							}
						}
					}
				}


					//System.out.println("topo : " + eElement.getElementsByTagName("Topology").item(0).getTextContent());

				}
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
		return oj;

	}

}