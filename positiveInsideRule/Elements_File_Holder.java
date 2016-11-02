package positiveInsideRule;

import java.util.ArrayList;

public class Elements_File_Holder {
	ArrayList <Elements_file> obj;
	public String file_name;
	int positive_score;
	int negative_score;
	
	Elements_File_Holder(){
		this.obj = new ArrayList <Elements_file> ();
		this.positive_score = 0;
		this.negative_score = 0;
	}

}
