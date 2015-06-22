package framework;

import java.util.Arrays;

/**
 * Stores the genomic interval of a TFBS
 * @author Thorsten Will
 */
public class BindingSite implements Comparable<BindingSite> {
	
	// stores the interval of the binding site
	private final int left;
	private final int right;
	
	// stores the orientation
	private final boolean plus_strand;

	BindingSite(int left, int right, boolean plus_strand) {
		this.left = left;
		this.right = right;
		this.plus_strand = plus_strand;
	}
	
	@Override
	public int compareTo(BindingSite bs) {
		
		int left_test = Integer.compare(this.left, bs.left);
		if (left_test == 0)
			return Integer.compare(this.right, bs.right);
		
		return left_test;
	}
	
	public String toString() {
		String sign = "+";
		if (!plus_strand)
			sign = "-";
		return "["+left+","+right+","+sign+"]";
	}
	
	/**
	 * Checks if BS are adjacent to each other according to given d_min/d_max
	 * returns 0 if incompatible (overlap/too far fro each other)
	 * returns -1 if other BS left of this BS -> [(l2,r2),(l1,r1)]
	 * returns +1 if this BS is left of other BS -> [(l1,r1),(l2,r2)]
	 * @param bs
	 * @return
	 */
	public int adjacentToRelation(BindingSite bs, int d_min, int d_max) {
		
		// this BS is right of the other BS
		if (this.left >= bs.left) {
			if ( this.left < (bs.right + d_min) || this.left > (bs.right + d_max) ) // overlap / too far away
				return 0;
			return -1;
		} else {
			if ( bs.left < (this.right + d_min) || bs.left > (this.right + d_max) ) // overlap / too far away
				return 0;
			return 1;
		}
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + left;
		result = prime * result + (plus_strand ? 1231 : 1237);
		result = prime * result + right;
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		BindingSite other = (BindingSite) obj;
		if (left != other.left)
			return false;
		if (plus_strand != other.plus_strand)
			return false;
		if (right != other.right)
			return false;
		return true;
	}

	public int getLeft() {
		return left;
	}

	public int getRight() {
		return right;
	}

	public boolean isPlusStrand() {
		return plus_strand;
	}
	
	/**
	 * Return readable String-repr. of an BS array
	 * @param bss
	 * @return
	 */
	public static String getBSArray(BindingSite[] bss) {
		return Arrays.asList(bss).toString();
	}
}
