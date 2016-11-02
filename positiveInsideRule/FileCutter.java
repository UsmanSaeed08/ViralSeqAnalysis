package positiveInsideRule;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;

import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojavax.SimpleNamespace;
import org.biojavax.bio.seq.RichSequence;
import org.biojavax.bio.seq.RichSequenceIterator;

public class FileCutter {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

		System.out.print("\n Cuts the fasta format protein input file to single sequences in the given directory \n");

		//String infile_path = "F:/yeast.fasta";
		//String outfile_path = "F:/yeastSplit";

		String infile_path = args[0];
		String outfile_path = args[1];
		if (infile_path.isEmpty() || outfile_path.isEmpty()){
			System.out.print("\n WrongInputParameters \n");
		}
		BufferedReader br ;
		BufferedWriter wr;

		try{
			br = new BufferedReader(new FileReader(infile_path)); //add file path here

			Alphabet alpha = AlphabetManager.alphabetForName("PROTEIN");
			SimpleNamespace ns = new SimpleNamespace("biojava");
			RichSequenceIterator iterator = RichSequence.IOTools.readFasta(br,
					alpha.getTokenization("token"), ns);

			int i = 0;
			while (iterator.hasNext()) {	//populate temp_file
				RichSequence rec = iterator.nextRichSequence();
				String acc = rec.getAccession();
				String seq = rec.seqString();
				String description = rec.getDescription();

				System.out.print(">"+acc+"| ");
				System.out.print(description+"\n");
				System.out.print(seq+"\n");
				i++;  
				wr = new BufferedWriter(new FileWriter(outfile_path+"_"+i+".fasta"));
				wr.write(">"+acc+"| ");
				if(!(description==null)){
					wr.write(description);
					wr.newLine();
				}
				else{
					wr.newLine();	
				}
				wr.write(seq);
				wr.newLine();

				//if (i>5){
					//break;
				//}
				wr.close();
			}
			br.close();

		}
		catch(Exception e){
			e.printStackTrace();
		}

	}

}
