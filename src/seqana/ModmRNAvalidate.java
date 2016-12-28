import acedb.*;
import java.io.*;
import java.util.regex.*;

public class ModmRNAvalidate {
	private Ace hgdb;

	/** right now AQL has a bug that prevented me from selecting the length of
	 *  the DNA sequence. 
	 * I have to iterate throught the whole database to get this number.
	 */
	private String query="select model::s, len, cDNA::s->corresponding_cDNA, s->corresponding_cDNA->DNA[2] from s in class sequence, len in (sum(select s->source_exons[2]) - sum (select s->source_exons[1]) + count(select s->source_exons)) where exists_tag s->corresponding_cDNA and len != s->corresponding_cDNA->DNA[2]";


	//private Result result;

	public ModmRNAvalidate() {
		hgdb = new Ace("gost", 3501, "kzhou", "fugufish");
	//	result = hgdb.execQuery(query);
	}
	public ModmRNAvalidate(String h, int p, String u, String pw) {
		hgdb = new Ace(h, p, u, pw);
	//	result = hgdb.execQuery(query);
	}

	public static void main(String[] args) {
		try {
			ModmRNAvalidate val = new ModmRNAvalidate();
			Result result = val.hgdb.execQuery(val.query);

			String mRNA, model;
			int mRNALen, modelLen, CDSend, cnt=0;
			//String seq;
			String fileSmall = "mRNAModSmall.diff";
			String filemRNALong = "mRNALong.diff";
			String filemRNAShort = "mRNAShort.diff";
			int small=0, mrnalong=0, mrnashort=0;

			PrintWriter outSmall = new PrintWriter( new BufferedWriter( new FileWriter(fileSmall)));
			PrintWriter outmRNAL = new PrintWriter( new BufferedWriter( new FileWriter(filemRNALong)));
			PrintWriter outmRNAS = new PrintWriter( new BufferedWriter( new FileWriter(filemRNAShort)));

			while (result.next()) {
				cnt++;
				if (cnt%100 == 0) 
					System.out.println("Working on ... #" + cnt + "\n" + result.getRowString());
				//out.println(result.getRowString());
				modelLen = result.getInt(2);
				mRNALen = result.getInt(4);
				if (Math.abs(modelLen-mRNALen) < 40) {
					outSmall.println(result.getRowString());
					small++;
					//out.println("small differences");
				}
				else if (mRNALen>modelLen) {
					Pattern pattern = Pattern.compile("a+\\z");
					// matching the polyA tail
					Matcher mm = pattern.matcher(val.hgdb.getDNA(result.getString(3))); 
					if (mm.find() && Math.abs(mRNALen - mm.group().length() - modelLen) < 40) {
						outSmall.print(result.getRowString());
						outSmall.println(" small differences by polyA: Ax" + mm.group().length());
						small++;
					}
					else {
						outmRNAL.println(result.getRowString());
						//outmRNAL.print(val.hgdb.getFastaDNA(result.getString(1))); // mRNA model 
						//outmRNAL.println(val.hgdb.getFastaDNA(result.getString(3)));  // cDNA (mRNA) 
						mrnalong++;
					}
				}
				else {
					outmRNAS.println(result.getRowString());
					//outmRNAS.print(val.hgdb.getFastaDNA(result.getString(1))); // mRNA model 
					//outmRNAS.println(val.hgdb.getFastaDNA(result.getString(3)));  // cDNA (mRNA) 
					mrnashort++;
				}
				/*
				model = result.getString(1);
				modelLen = result.getInt(2);
				mRNA = result.getString(3);
				mRNALen = result.getInt(4);
				if (modelLen != mRNALen) {
					System.out.println(model + ": " + modelLen + " != " + mRNA + ": " + mRNALen);
				}
				*/
			}
			outSmall.close();
			outmRNAL.close();
			outmRNAS.close();
			System.out.println("small diff: " + small + " mRNA Long: " + mrnalong + " mRNA Short: " + mrnashort);
		}
		catch (AceException ae) {
			System.err.println(ae);
			System.exit(1);
		}
		catch(Exception e) {
			System.err.println(e);
			System.exit(1);
		}
	}
}
