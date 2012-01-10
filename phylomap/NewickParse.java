package phylomap;
import java.io.*;
import java.util.*;

public class NewickParse {

	/**
	 * @param args
	 */
	
	public node parseFromFile(File fin) throws Exception{
		node root=null;
		Stack stack=new Stack();
		FileReader fr=new FileReader(fin);
		BufferedReader bf=new BufferedReader(fr);
		String s=bf.readLine();
		String tree="";
		while(s!=null){
			tree+=s.trim();
			s=bf.readLine();
		}
		
		bf.close();
		fr.close();
		
		// store node_name : node reference 
		HashMap hm=new HashMap();
		
		String nextChar=tree.substring(0, 1);
		int charCounter=0;
		int p_node_counter=0;
		while(nextChar!=null && nextChar!=";"){
			
			if(")".equals(nextChar)){
				String nodeString="";
				String nextPop=(String)stack.pop();
				while(!nextPop.equals("(")){
					nodeString=nextPop+nodeString;
					nextPop=(String)stack.pop();
				}
				//parse node name and length
				String[] nodes=nodeString.split(",");
				if(nodes.length==2){
					String nodeA=nodes[0];
					String nodeB=nodes[1];
					String[] nodeAs=nodeA.split(":");
					String[] nodeBs=nodeB.split(":");
					String nodeA_name=nodeAs[0];
					String nodeA_value=nodeAs[1];
					String nodeB_name=nodeBs[0];
					String nodeB_value=nodeBs[1];
					
					//check for if node already exist, otherwise create a new node for parsed one
					node Node_A=new node();
					double dnode_a=Double.parseDouble(nodeA_value);
					if(nodeA_name.startsWith("#")){
						Node_A=(node)hm.get(nodeA_name);
						Node_A.setDistance(dnode_a);
					}else{
						Node_A.setId(nodeA_name);
						Node_A.setDistance(dnode_a);
					}
					
					
					node Node_B=new node();
					double dnode_b=Double.parseDouble(nodeB_value);
					if(nodeB_name.startsWith("#")){
						Node_B=(node)hm.get(nodeB_name);
						Node_B.setDistance(dnode_b);
					}else{
						Node_B.setId(nodeB_name);
						Node_B.setDistance(dnode_b);
					}
					
					//create a parent node to connect both
					node Node_P=new node();
					Node_P.setId("#"+p_node_counter++);
					Node_P.setLeft(Node_A);
					Node_P.setRight(Node_B);
					hm.put(Node_P.getId(), Node_P);
					
					// push the new node name into stack
					for(int i=0;i<Node_P.getId().length();i++){
						stack.push(Node_P.getId().substring(i, i+1));
					}
				}else{
					String nodeA=nodes[0];
					String nodeB=nodes[1];
					String nodeC=nodes[2];
					String[] nodeAs=nodeA.split(":");
					String[] nodeBs=nodeB.split(":");
					String[] nodeCs=nodeC.split(":");
					String nodeA_name=nodeAs[0];
					String nodeA_value=nodeAs[1];
					String nodeB_name=nodeBs[0];
					String nodeB_value=nodeBs[1];
					String nodeC_name=nodeCs[0];
					String nodeC_value=nodeCs[1];
					
					//check for if node already exist, otherwise create a new node for parsed one
					node Node_A=new node();
					double dnode_a=Double.parseDouble(nodeA_value);
					if(nodeA_name.startsWith("#")){
						Node_A=(node)hm.get(nodeA_name);
						Node_A.setDistance(dnode_a);
					}else{
						Node_A.setId(nodeA_name);
						Node_A.setDistance(dnode_a);
					}
					
					
					node Node_B=new node();
					double dnode_b=Double.parseDouble(nodeB_value);
					if(nodeB_name.startsWith("#")){
						Node_B=(node)hm.get(nodeB_name);
						Node_B.setDistance(dnode_b);
					}else{
						Node_B.setId(nodeB_name);
						Node_B.setDistance(dnode_b);
					}
					
					node Node_C=new node();
					double dnode_c=Double.parseDouble(nodeC_value);
					if(nodeC_name.startsWith("#")){
						Node_C=(node)hm.get(nodeC_name);
						Node_C.setDistance(dnode_c);
					}else{
						Node_C.setId(nodeC_name);
						Node_C.setDistance(dnode_c);
					}
					
					//create a parent node to connect both
					node Node_P=new node();
					Node_P.setId("#"+p_node_counter++);
					Node_P.setLeft(Node_A);
					Node_P.setRight(Node_B);
					Node_P.setDistance(dnode_c);
					hm.put(Node_P.getId(), Node_P);
					
					node Node_root=new node();
					Node_root.setId(nodeC_name);
					Node_root.setLeft(Node_P);
					Node_root.setRight(null);
					//Node_root.setDistance(dnode_c);
					Node_root.setInfo("root");
					hm.put(Node_root.getId(), Node_root);
					
					// push the new node name into stack
					for(int i=0;i<Node_root.getId().length();i++){
						stack.push(Node_root.getId().substring(i, i+1));
						
					}
				}
			}else{
				stack.push(nextChar);
				
			}
			charCounter++;
			if(charCounter>tree.length()-1){
				nextChar=null;
			}else{
				nextChar=tree.substring(charCounter, charCounter+1);
			}
		}
		
		String root_s="";
		stack.pop();
		while(!stack.empty()){
			root_s=(String)stack.pop()+root_s;
		}
		root=(node)hm.get(root_s);
		
		return root;
	}
	
	
	public static void main(String[] args) {
		
		File treefile=new File("/storage/disk2/NS1.tree");
		File fout=new File("F:\\matlab\\NG\\tree_dis.dat");
		NewickParse parser = new NewickParse(); 
		
		try{
			node root=parser.parseFromFile(treefile);
			postOrder po=new postOrder();
			System.out.println("Post order the tree................................");
			System.out.println("root:"+root.getId()+":"+root.getDistance());
			//po.genDisMatrixFromTree(root, fout);
			po.postTree(root.getLeft());
			//po.postTree_cal_tree_dis(root.getLeft());
			//po.printHm();
			//System.out.println(root.getLeft().getId()+":"+root.getLeft().getDistance());
			
		}catch(Exception e){
			e.printStackTrace();
		}
	}

}
