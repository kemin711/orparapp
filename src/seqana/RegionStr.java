public class RegionStr extends Region 
{
	// simply add a string for identifier for this region
	// use RegionInt if you want to use a integer to identify a
	// region.

	private String key;

	RegionStr(String name, int s, int e) {
		super(s,e);
		key=name;
	}
}
