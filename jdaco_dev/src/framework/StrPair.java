package framework;

/**
 * Helper class, stores sorted par of strings
 * @author Thorsten Will
 */
public class StrPair {
	private final String l;
	private final String r;
	
	public StrPair(String l, String r) {
		if (l.compareTo(r) < 0) {
			this.l = l;
			this.r = r;
		} else {
			this.l = r;
			this.r = l;
		}
	}

	@Override
	public int hashCode() {
		return l.hashCode() + r.hashCode();
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		StrPair other = (StrPair) obj;
		if (l == null) {
			if (other.l != null)
				return false;
		} else if (!l.equals(other.l))
			return false;
		if (r == null) {
			if (other.r != null)
				return false;
		} else if (!r.equals(other.r))
			return false;
		return true;
	}

	public String getL() {
		return l;
	}

	public String getR() {
		return r;
	}

	@Override
	public String toString() {
		return l +" "+ r;
	}
	
}
