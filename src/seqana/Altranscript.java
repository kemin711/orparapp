import acedb.*;
import java.io.*;

public class Altranscript 
{
	public static void fileinput(String inf) {
		try {
			FileReader fileread = new FileReader(inf);
			BufferedReader br = new BufferedReader(fileread);

			String line;
			String[] words;
			line = br.readLine();
			int genomicLen, start, end;
			String genomicName, subName;
			if (line.charAt(0) == '"') line = line.substring(1);
			words = line.split("\"?\t\"?");
			RegionStr reg;
			while (line != null) {
				System.out.println();
				genomicName = words[0];
				genomicLen = Integer.parseInt(words[1]);
				System.out.println(genomicName + " len=" + genomicLen);
				while (line != null && words[0].equals(genomicName) && 
						genomicLen == Integer.parseInt(words[1]) )
				{
					System.out.println(words[2] + " || " + words[3] + " || " + words[4]);
					reg = new RegionStr(words[2], Integer.parseInt(words[3]), Integer.parseInt(words[4]));

					line = br.readLine();
					if (line != null) {
						if (line.charAt(0) == '"') line = line.substring(1);
						words = line.split("\"?\t\"?");
					}
					else break;
				}
			}
			br.close();
			fileread.close();
		}
		catch (FileNotFoundException fe) {
			System.err.println("File not found " + fe);
			System.exit(1);
		}
		catch (IOException e) {
			System.err.println(e);
			System.exit(1);
		}
	}
	/** input from acedb with aql query has a bug in acedb.
	 *  it cannot be donw right now. I have to use tableMaker
	 *  to extract the table first into a file.  Then read from
	 *  file.
	 */
	public static void dbinput() {
		try {
			Ace hsg = new Ace("gost", 3501, "kzhou", "fugufish");
			//String aql = "select s from s in class sequence where exists_tag s->Genomic_canonical";
			String aql = "select s, s->DNA[2] from s in class sequence where exists_tag s->Genomic_canonical";

			Result rslt = hsg.execQuery(aql);
			String genid;
			int count = 0;
			while (rslt.next()) {
				count++;
				genid = rslt.getString(1);
				int seqlen = rslt.getInt(2);
				System.out.print("\n" + count + "=>" + genid + "  " + seqlen + "\n");
				//System.out.print(genid + "\n");
				//String subaql = "select s, ss->DNA[2] from s in object(\"Sequence\", \"" + genid +"\")->Subsequence, ss in class Sequence where ss = \"" + genid + "\"";
				//String subaql = "select s from s in object(\"Sequence\", \"" + genid +"\")->Subsequence";
				String subaql = "select s->Subsequence, s->Subsequence[2], s->Subsequence[3] from s in class Sequence where s = \"" + genid +"\"";
				//Result rsltSub = hsg.execQuery("select s, s->DNA[2] from s in object(\"Sequence\", " + genid +")->subsequence");
				try {
					Result rsltSub = hsg.execQuery(subaql);
					if (rsltSub.isEmpty()) {
						System.err.println(genid + " has no Subsequences");
					}
					else {
						while (rsltSub.next()) {
							System.out.print(rsltSub.getRowString() + "\n");
						}
					}
				}
				catch (AceException sae) {
					System.err.println(sae);
					System.err.println(subaql);
					System.exit(1);
				}
			}
		}
		catch (AceException ae) {
			ae.printStackTrace();
			System.exit(1);
		}
	}
	public static void main(String[] args) {
		//dbinput();
		fileinput("/data/genome/hs/genSubseq.txt");
	}
}
