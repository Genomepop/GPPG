// Configuration options
{
	"individuals": 10000,
	"generations": 10000,
	"scaling": 100,
	"compression": {"name": "Greedy-Load", "k": 20, "t":10},
	//"compression": {"name": "Store-Active", "k": 10},
	"genotype": {"name": "Sequence", "length":10000	},
	"operators": [
				  {"name":"PointMutation", "rate":1e-9},
				  {"name":"Insertion", "rate":1e-9, "size":[10,20]},
				  {"name":"Deletion", "rate":1e-9, "size":[10,20]},
				  {"name":"Recombination", "rate":1e-9}
				  ],
	//"genotype": {"name": "Pathway", "genes":100, "tfs":25, "regions":[10,100]}
	//{"name":"BindingSiteMutation", "overlap":5, "gainRate":1e-10, "lossRate":[1e-9,0.7] }
	"output": {"individuals":false, "performance":"path/to/file"}
}