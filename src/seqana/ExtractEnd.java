import acedb.*;
import java.util.TreeMap;
import java.util.Vector;
import java.io.*;

/** ExtractEnd class 
 * extract the 3'-end of those sequences whose
 *	3'-end of CDS and mRNA are the same.
 * each object should be based on one genome sequence.  The main method
 * iterate through the main query result.
 * The genomic sequence is too long in this case. Combined with the 
 * multiple threads running and the poor java efficiency (The unocode
 * is one factor.  The subsequence extraction cannot be done from this
 * computer.  The output of this program is the instruction files as to
 * how to obtain the subsequences.  
 * Then I need to write a separate perl program (exend.pl)
 *  to extract the subsequences.
 */

public class ExtractEnd {
	private String gseq;        // name of the genomic sequence
//	private String DNA;         // the genomic DNA, at testing stage, no real extraction yet
	private int DNALength;
	private TreeMap name2regv;   // name->locsorted index
	private Vector locsorted;
	private Ace ace;           // human genome acedb
	private PrintWriter out;

	public final int TAIL_LIMIT = 5000;    // 5 kb 3'-end max extracted

	/** basic set up, not ready for action yet.
	  */
	public ExtractEnd(Ace db, String outf) {
		ace = db;
		name2regv = new TreeMap();
		locsorted = new Vector(100);
		try {
			out = new PrintWriter(new BufferedWriter(new FileWriter(outf)));
		}
		catch (IOException e) {
			System.err.println(e);
			System.exit(1);
		}
	}

	/** set up this object for extraction.
	  * It first clears the internal map and vector member.
	  */
	public void setup(String gn, int len) {
		gseq = gn;
		DNALength = len;
		name2regv.clear();  // this is more efficient, no reallocation needed
		locsorted.clear();  // no reallocation of memory

		try {
			Result subrslt = ace.execQuery("select subseq::s, start::s[1], end::s[2] from s in object(\"Sequence\",\"" + gseq + "\")->Subsequence order by :2, :3");

			//System.out.println(gseq);
			int i=0;  
			while (subrslt.next()) {
				locsorted.add(new Region(subrslt.getInt(2), subrslt.getInt(3)));
				name2regv.put(subrslt.getString(1), new Integer(i++));
			}
//			System.err.println("before getting DNA sequence");
//			DNA = ace.getDNA(gseq);  // at debug stage, this one is not called yet
			// running out of memory on this computer.
//			System.err.println("after getting DNA sequence");
		}
		catch (AceException ae) {
			System.err.println(ae);
			System.exit(1);
		}
		catch (Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
	}

	public void extract() {
		try {
			Result subrslt = ace.execQuery("select s, s->CDSofmodel from s in object(\"Sequence\", \"" + gseq + "\")->subsequence where exists_tag s->cDNA_model and exists s->CDSofmodel");
			String cds, model;
			int modelEnd, cdsEnd;
			int i;    // the index into the locsorted vector

			Region modelR, cdsR, nextR, nextOverlapR;

			//out.println("Genomic " + gseq);  // may not have anything to show
			// delayed untill knowing the size of the buffer
			StringBuffer buff = new StringBuffer(1000);
			while (subrslt.next()) {  // used the model; cds would yield the same result
				model = subrslt.getString(1);
				cds = subrslt.getString(2);
				i = ((Integer)name2regv.get(model)).intValue();
				modelR = (Region)locsorted.get(i);
				cdsR = (Region)locsorted.get( ((Integer)name2regv.get(cds)).intValue() );

				if (modelR.differentEnd(cdsR)) continue; 

				modelEnd = modelR.getEnd();
				//out.print(">" + model + " 3'-region\t");
				buff.append(">" + model + " 3'-region\t");
				if (modelR.forward()) {
               i++;
               while (true) {
                  while (i<locsorted.size() && modelR.overlap((Region)locsorted.get(i))>1 )
                     i++;
                  if (i == locsorted.size()) { // reached the end of the genomic sequence
                     //if (DNA.length()-modelEnd > TAIL_LIMIT)
                     if (DNALength-modelEnd > TAIL_LIMIT)
                        //Ace.writeFasta(out, DNA.substring(modelEnd, modelEnd+TAIL_LIMIT), 70);
                        //out.println(modelEnd + " " + TAIL_LIMIT);  // startIndex, length
                        buff.append(modelEnd + " " + TAIL_LIMIT + "\n");  // startIndex, length
                     //else Ace.writeFasta(out, DNA.substring(modelEnd), 70);
                     else buff.append(modelEnd + "\n");   // all the way to the end
                     //else // out.println(modelEnd);   // all the way to the end
                     //System.out.println(model + " at the end of genomic sequence\n");
                     break;
                  }
                  else {
                     nextR = (Region)locsorted.get(i);
                     i++;
                     while (i<locsorted.size() && nextR.overlap(nextOverlapR=(Region)locsorted.get(i))>0 ) {
                        if (modelR.distance(nextOverlapR) < modelR.distance(nextR))
                           nextR = nextOverlapR;
                        i++;
                     }
                     if (modelR.overlap(nextR) > 0) continue;
                     else {
                        //Ace.writeFasta(out, DNA.substring(modelEnd, Math.min(modelEnd+TAIL_LIMIT, nextR.getSmallEnd())), 70);
                        //out.println(modelEnd + " " + Math.min(TAIL_LIMIT, modelR.distance(nextR)));
                        buff.append(modelEnd + " " + Math.min(TAIL_LIMIT, modelR.distance(nextR)) + "\n");
                        // forward direction
                        //System.out.println("=>get subseq " + model + " " + modelR + " to " + nextR);
                        break;
                     }
                  }
               } // while loop
            }
            else { // gene backward direction
               i--;
               while (true) {
                  while (i>=0 && modelR.overlap((Region)locsorted.get(i))>1 )
                     i--;
                  if (i<0) { // first gene on the genomic sequence, reached the beginning
                     //Ace.writeFasta(out, Ace.reverseComplement(DNA.substring(Math.max(modelEnd-TAIL_LIMIT-1, 0), modelEnd-1)) ,70);
                     //if (modelEnd-TAIL_LIMIT-1 <= 0) out.print("0 " + (modelEnd-1));
                     if (modelEnd-TAIL_LIMIT-1 <= 0) buff.append("0 " + (modelEnd-1));
                     else buff.append(modelEnd-TAIL_LIMIT-1 + " " + TAIL_LIMIT);
                     //else out.print(modelEnd-TAIL_LIMIT-1 + " " + TAIL_LIMIT);
                     //System.out.println(model + " is the firt gene on the genomic DNA");
                     buff.append(" rc\n");
                     break;
                  }
                  else {
                     nextR = (Region)locsorted.get(i);

                     i--;
                     while (i>=0 && nextR.overlap(nextOverlapR=(Region)locsorted.get(i))>0 )
                     {
                        if (modelR.distance(nextOverlapR) < modelR.distance(nextR))
                           nextR = nextOverlapR;
                        i--;
                     }
                     if (modelR.overlap(nextR)>0) continue;
                     else {
                        //backward direction
                        //Ace.writeFasta(out, Ace.reverseComplement( DNA.substring(Math.max(nextR.getLargeEnd(), modelEnd-TAIL_LIMIT-1), modelEnd)) ,70);
                        if (modelR.distance(nextR) <= TAIL_LIMIT)
                           buff.append(nextR.getLargeEnd() + " " + (modelR.distance(nextR)-1));
                           //out.print(nextR.getLargeEnd() + " " + (modelR.distance(nextR)-1));
                        else buff.append(modelEnd-TAIL_LIMIT-1 + " " + TAIL_LIMIT);
                        //else out.print(modelEnd-TAIL_LIMIT-1 + " " + TAIL_LIMIT);
                        //System.out.println("<=get subseq " + model + " " + modelR + " to " + nextR);
                        buff.append(" rc\n");
                        break;
                     }
                  }
                  //out.println(" rc");
                  //buff.append(" rc\n");
               } // while loop
            } // backward direction
         }
         if (buff.length()>0) out.print("Genomic " + gseq + "\n" + buff.toString());
		}
		catch (AceException ae) {
			System.err.println(ae);
			System.exit(1);
		}
	}

	public void close() {
		out.close();
		ace.close();
	}

		
	public static void main(String[] args) {
		try {
			Ace hsgn = new Ace("gost", 3501, "kzhou", "fugufish");
			String query="select s, s->DNA[2] from s in class sequence where exists_tag s->genomic_canonical and exists_tag s->subsequence";
			Result genomicResult = hsgn.execQuery(query);
			ExtractEnd exnd = new ExtractEnd(hsgn, "tail.ext");

			int count=0;
			while (genomicResult.next()) {
				System.out.println(++count + " " + genomicResult.getString(1));
				exnd.setup(genomicResult.getString(1), genomicResult.getInt(2));
				exnd.extract();
			}
			System.out.println("result written to tail.ext file. all done!");
			exnd.close();
		}
		catch (AceException ae) {
			System.err.println(ae);
			ae.printStackTrace();
			System.exit(1);
		}
		catch (Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
	}
}
