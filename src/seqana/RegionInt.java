public class RegionInt extends Region 
{
	// simply add a string for identifier for this region
	// use RegionInt if you want to use a integer to identify a
	// region.
	// or I have to use object as the key and do nasty casting

	private int key;

	RegionInt(int name, int s, int e) {
		super(s,e);
		key=name;
	}
}
