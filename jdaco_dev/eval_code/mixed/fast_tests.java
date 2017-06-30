package mixed;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import framework.Utilities;

public class fast_tests {
	
	public static void main(String[] args) {
		Map<String, List<String>> graph = new HashMap<>();
		
		List<String> a_list = new ArrayList<>();
		a_list.add("b");
		a_list.add("c");
		graph.put("x", a_list);
		
		List<String> b_list = new ArrayList<>();
		b_list.add("x");
		graph.put("b", b_list);
		
		List<String> c_list = new ArrayList<>();
		c_list.add("d");
		c_list.add("e");
		c_list.add("f");
		graph.put("c", c_list);
		
		List<String> d_list = new ArrayList<>();
		d_list.add("b");
		graph.put("d", d_list);
		
		List<String> e_list = new ArrayList<>();
		e_list.add("h");
		e_list.add("g");
		graph.put("e", e_list);
		
		List<String> f_list = new ArrayList<>();
		f_list.add("g");
		graph.put("f", f_list);
		
		System.out.println(graph);
		System.out.println(Utilities.getSCCs(graph));
	}
}
