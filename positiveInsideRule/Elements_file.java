package positiveInsideRule;

import java.util.ArrayList;

// ds for handling the file elements and sequences etc
public class Elements_file {
	
	public String id;
	public int tm_no;
	public boolean SP;
	int sp_begin;
	int sp_end;
	public String file_name;
	public String seq;
	public ArrayList<String> tm_segments;
	public int[] start ;
	public int[] end = new int [10];
	public float[] segment_scorep;		//initialized with tm_no and index->tmSegment & value->Score
	public float[] segment_scoren;
	public int scorep;	//gaps in tm region - total accumulation from all 8 hits
	public int scoren;	//positive residue mismatches in tm - total accumulation from all 8 hits 
	public float coverage_querytm;
	public float coverage_Hittm;
	public String md5;
	public int nterm; //1-inside and 0-outside
	public int cterm; //1-inside and 0-outside
	
	public Elements_file(){
		this.nterm = 3;
		this.cterm = 3;
		this.SP = false;
		this.sp_begin = 0;
		this.sp_end = 0;
		this.coverage_querytm = 0f;
		this.coverage_Hittm = 0f;
		this.id = "";
		this.tm_no = 0;
		this.start = new int [200];
		this.end = new int [200];
		this.tm_segments = new ArrayList <String> ();
		this.seq = "";
		this.scorep = 0;
		this.scoren = 0;
		this.md5 = "";
	}
	
	public void intializeScore (int v){
		this.segment_scorep = new float[v];
		this.segment_scoren = new float[v];
	}

}
