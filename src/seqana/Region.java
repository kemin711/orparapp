//package Region;

public class Region implements Comparable
{
	private int start;
	private int end;
	//private String key;  // a name for the region, should be unique
	// this limits the key to be a string, in reality, a key can be
	// an integer. So I am not sure whether this is good or not.

	/*
	static final int FORWARD=1;
	static final int BACKWARD=-1;
	static final int NODIRECTION = 0;
	*/

	/** constructor, create a region from s (start) and e (end).
	 */
	public Region(int s, int e) {
		start = s;
		end =e;
	}

	/** return the directionality of the region.
	 * @return end-start as an integer.
	 */
	public int direction() { 
		return end-start; 
	}
	public boolean backward() { return end<start; }
	public boolean forward() { return end>start; }

	public String toString() {
		return new String(start + "," + end);
	}

	/** return the smaller end regardless of direction.
	 */
	public int getSmallEnd() {
		if (end>start) return start;
		else return end;
	}
	public int getLargeEnd() {
		if (end<start) return start;
		else return end;
	}
	public int getEnd() { return end; }

	/** return the distance of this Region to r.
	 * Returns a negative number if overlap.
	 */
	public int distance(Region r) {
		// --this-->
		if (end>start) {
			if (r.end>r.start) { // both forward
				//-this-> -r-> or --r--> --this-->
			 	if (end < r.start) return r.start - end;
			 	else if (r.end < start) return start - r.end;
			  	else return Math.max(start,r.start)-Math.min(end,r.end);
			}
			// <--r--
			else {
				if (end<r.end) return r.end-end;
				else if (start > r.start) return start-r.start;
				else return Math.max(start, r.end)-Math.min(end, r.start);
			}
		}
		else {
			if (r.end>r.start) { //<-this- -r->
				if (start < r.start) return r.start-start;
				else if (end > r.end) return end-r.end;
				else return Math.max(end,r.start)-Math.min(start,r.end);
			}
			else {
				if (start < r.end) return r.end- start;
				else if (r.start < end) return end - r.start;
				else return Math.max(end,r.end)-Math.min(start,r.start);
			}
		}
	}

	public int compareEnd(Region r) {
		return end - r.end;
	}
	public boolean differentEnd(Region r) {
		return end != r.end;
	}
	public boolean sameEnd(Region r) {
		return end==r.end;
	}

	/** return true if this region and r have the same start and end.
	 */
	public boolean equals(Region r) {
		return start==r.start && end==r.end;
	}
	public int hashCode() { // for hashMap
		return start + end%113;  // at most 100 alternative transcript ends
	}

	public int compareTo(Object r) {  // for sorting
		int result = ((Region)r).start-start;
		return result==0? ((Region)r).end-end : result;
	}

	/** return the extent of overlap between this Region and r.
	 * In case you want to know the extent of the overlap, this
	 * method is very useful.
	 * @return if overlap returns a positive number.
	 *         a negative number indicates no overlap.
	 */
	public int directionalOverlap(Region r) {
		if (end > start && r.end > r.start) {
			// s-->e
			//       r.s-->r.e
			if (r.start > end) return end-r.start+1;
			//              s-->e
			// r.s-->r.e
			else if (start > r.end) return r.end-start+1;
			else {
				// s----->e      s------->e    s-->e         s----->e
				//  r.s---->r.e   r.s-->r.e r.s----->r.e r.s---->r.e
				return Math.min(end,r.end) - Math.max(start,r.start) + 1;
			}
		}
		else if (end < start && r.end < r.start) {
			// e<---s 
			//        r.e<---r.s
			if (r.end > start) return start-r.end+1;
			//              e<---s 
			// r.e<---r.s
			else if (end > r.start) return r.start-end+1;
			//              e<---s 
			//          r.e<---r.s
			else return Math.max(start,r.start)-Math.min(end,r.end)+1;
		}
		else return -1;
	}

	/** return the overlap length.
	 * @return a negative value indicates no overlap;
	           a positive value indicates a overlap.
	*/
	public int overlap(Region r) 
	{
		if (end > start && r.end > r.start) { // both forward
			// s-->e
			//       r.s-->r.e
			if (r.start > end) return end-r.start+1;
			//              s-->e
			// r.s-->r.e
			else if (start > r.end) return r.end-start+1;
			else {
				// s----->e      s------->e    s-->e         s----->e
				//  r.s---->r.e   r.s-->r.e r.s----->r.e r.s---->r.e
				return Math.min(end,r.end) - Math.max(start,r.start) + 1;
			}
		}
		else if (end > start  && r.end < r.start) { // forward, backward
			// s------>e
			//          r.e<------r.s
			if (r.end > end) return end-r.end+1;
			//                s------>e
			// r.e<------r.s
			else if (start > r.start) return r.start-start+1;
			// s------>e
			//    r.e<------r.s
			else return Math.max(end,r.start) - Math.min(start,r.end)+1;
		}
		else if (end < start && r.end > r.start) { //backward,forward
			// e<---s
			//        r.s--->r.e
			if (start < r.start) return start-r.start+1;
			//            e<---s
			// r.s--->r.e
			else if (end > r.end) return r.end-end+1;
			//            e<---s
			//        r.s--->r.e
			else return Math.max(start,r.end)-Math.min(end,r.start)+1;
		}
		else if (end < start && r.end < r.start) { // both backward
			// e<---s 
			//        r.e<---r.s
			if (r.end > start) return start-r.end+1;
			//              e<---s 
			// r.e<---r.s
			else if (end > r.start) return r.start-end+1;
			//              e<---s 
			//          r.e<---r.s
			else return Math.max(start,r.start)-Math.min(end,r.end)+1;
		}
		else { // dir == 0 || r.dir == 0
			// at most one base overlap
			if (start == end){
				System.out.println("it is a point: " + this);
				System.exit(1);
				if (start >= r.start && start <= r.end) return 1;
				else if (start < Math.min(r.start,r.end)) 
					return start - Math.min(r.start,r.end) + 1;
				else // (start > Math.max(r.start,r.end))
					return start - Math.max(r.start,r.end) + 1;
			}
			else {
				// (r.start == r.end) {
				System.out.println("this: " + this + " compared to target is a point: " + r);
				System.exit(1);
				if (r.start >= start && r.start <= end) return 1;
				else if (r.start < Math.min(start,end)) 
					return r.start-Math.min(start,end) + 1;
				else //if (r.start > Math.max(start,end))
					return r.start-Math.max(start,end)+1;
			}
		}
	}

	public static void main(String[] args) 
	{
		Region subseq = new Region(300, 5000);
		Region subseq2 = new Region(500, 6000);
		System.out.print("Testing " + subseq + " and " + subseq2);
		System.out.print("\nOverlap is " + subseq.overlap(subseq2));
		System.out.print("\nend\n");
	}
}
