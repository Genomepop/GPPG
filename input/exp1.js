// Configuration options
{
	"individuals": 1000,
	"generations": 1000,
	"scaling": 1,
	"steps":100,
	"compression": {"name": "Greedy-Load", "k": 50, "t":5},
	//"compression": {"name": "Store-Active"},
	"genotype": {"name": "Sequence", "length":100000 },
	"operators": [
				  {"name":"PointMutation", "rate":1e-4, "cost":1},
				  {"name":"Insertion", "rate":1e-4, "cost":1, "size":[10,100]},
				  {"name":"Deletion", "rate":1e-4, "cost":1, "size":[10,100]}
				  ],
	"output": {"individuals":"path/to/folder", "performance":"/Users/truths/tmp/exp1.csv"}
}