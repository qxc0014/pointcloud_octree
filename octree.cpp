#include<iostream>
#include<fstream>
#include<vector>
#include<pcl-1.9/pcl/io/pcd_io.h>
#include<pcl-1.9/pcl/filters/voxel_grid.h>
#include<pcl-1.9/pcl/point_types.h>
#include<pcl-1.9/pcl/visualization/cloud_viewer.h>
#include<pcl-1.9/pcl/common/transforms.h>
#include <chrono>
using namespace std;
using namespace pcl;
const float min_extent = 1;
const float com_float = 2e6;
int leafsize =0; 
string db_list = "/home/esoman/c++code/c++/octree/000000.bin";
class Octant{
    public:
        Octant(){}
        Octant(Octant* child_root,PointXYZ center,float extent,vector<int> points_index,bool is_leaf,int depth)
        :extent_(extent),points_index_(points_index),is_leaf_(is_leaf),center_(center),depth_(depth){
           // for(int i=0;i<8;i++)
           //    child_root_[i] = new Octant();
        }
        int depth_=0;//节点的深度
        Octant* child_root_[8]={nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr};//存八个子立方体的指针
        PointXYZ center_;//当前立方体的中心坐标
        float extent_;//当前立方体的半边长
        vector<int> points_index_;//当前立方体的包含点的Index
        bool is_leaf_;//当前坐标是否为叶子
        
};
class distindex{
    public:
        distindex(float dist_,int index_):dist(dist_),index(index_){}
        float dist;
        int index;
};
class result{
    public:
    result(float worst_dis_):worst_dis(worst_dis_){}//用于搜索一个近邻点
    result(float worst_dis_,int k):worst_dis(worst_dis_),worst_dis_cap(vector<distindex>(k,distindex(worst_dis_,-1))),size(k){ 
    }
    float worst_dis=0;
    int index;
    int num=0;
    int size;
    vector<distindex> worst_dis_cap;
    void add_point(float bias,int node_index);
};
void result::add_point(float bias,int node_index){
        if(num != this->size) num++;//已插入值的个数
        if(bias >= worst_dis_cap[this->size-1].dist) return;//大于最大值直接跳出
        int i = num-1;//已经插入最大值的index
        while(i>0){
             if(bias < worst_dis_cap[i-1].dist){
                this->worst_dis_cap[i] = worst_dis_cap[i-1];
                i--;
             }else{
                break;
             }
        }
        worst_dis_cap[i].dist = bias;
        worst_dis_cap[i].index = node_index;
        this->worst_dis =  worst_dis_cap[this->size-1].dist;
}
Octant* build_octree(Octant* root,pcl::PointCloud<PointXYZ>::Ptr db,PointXYZ center,float extent,vector<int> points_index,int depth,int width){
    if(points_index.size() == 0) {
        return nullptr;
    }
    if(root == nullptr){
        depth++;
       // cout << "节点深度：" << depth << "节点宽度" << width << endl;
        root = new Octant(nullptr,center,extent,points_index,true,depth);
    }
    if(extent < min_extent && points_index.size()<=1){
        root->is_leaf_ = true;//叶子节点
    }else{
        root->is_leaf_ = false;//不是叶子
        vector<vector<int>> child_point_index(8);
        for(auto point_idx:points_index){
            int Coordinate = 0;
            if((*db)[point_idx].x > center.x){
                Coordinate = Coordinate | 1;
            }
            if((*db)[point_idx].y > center.y){
                Coordinate = Coordinate | 2;
            }
            if((*db)[point_idx].z > center.z){
                Coordinate = Coordinate | 4;
            }
            child_point_index[Coordinate].push_back(point_idx);
        }
        float factor[2] = {-0.5,0.5};
        vector<PointXYZ> child_center(8);
        float child_extent=0;
        for(int i = 0;i < 8;i++){
            child_center[i].x = center.x + factor[(i&1)>0]*extent;
            child_center[i].y = center.y + factor[(i&2)>0]*extent;
            child_center[i].z = center.z + factor[(i&4)>0]*extent;
            child_extent = 0.5 *extent;
            //cout << child_extent << endl;
            root->child_root_[i] = build_octree(root->child_root_[i],db,child_center[i],child_extent,child_point_index[i],depth,i);
        }   
    }
    return root;
} 
//判断球与立方体的方位
bool overlap(Octant* root,PointXYZ Point,float worst_dis){
    //分三种情况:
    //第一种:球与立方体没有接触,只要投影的某个方向满足就可以
    float xyz[3];
    xyz[0] = fabs(root->center_.x - Point.x);
    xyz[1] = fabs(root->center_.y - Point.y);
    xyz[2] = fabs(root->center_.z - Point.z);
    float max_dis = (root->extent_+ worst_dis);
    if( xyz[0] > max_dis || xyz[1] > max_dis || xyz[2] > max_dis) return false;
    //第二种:球与立方体相交（通过投影判断）至少有两个投影面包含了圆心就可以认为是相交
    if(((xyz[0]<root->extent_)+(xyz[1]<root->extent_)+(xyz[2]<root->extent_))>=2) return true;
    //第三种:补充第二种，在边界处相交不满足第二种
    float x = (xyz[0]-root->extent_)>0?(xyz[0]-root->extent_):0;
    float y = (xyz[1]-root->extent_)>0?(xyz[1]-root->extent_):0;
    float z = (xyz[2]-root->extent_)>0?(xyz[2]-root->extent_):0;
    if(x*x+y*y+z*z<worst_dis*worst_dis) return true;
}
//判断球是否在立方体内
bool inside(Octant* root,PointXYZ Point,float worst_dis){
    float xyz[3];
    xyz[0] = fabs(root->center_.x - Point.x);
    xyz[1] = fabs(root->center_.y - Point.y);
    xyz[2] = fabs(root->center_.z - Point.z);
    float max_dis = (root->extent_ - worst_dis);
    return ((xyz[0] < max_dis) && (xyz[1] < max_dis) && (xyz[2] < max_dis));

}
bool octree_knn_search(Octant* root,pcl::PointCloud<PointXYZ>::Ptr db,PointXYZ Point,result &a){
    //先判断当前root是否为空指针
    if(root == nullptr) return false;
    //判断当前的节点是否为叶子
    if((root->is_leaf_ == true) && root->points_index_.size() == 1){
       //计算worst_dis
       //cout << "找到叶子！" << endl;
       Eigen::Vector3d radius(Point.x - (*db)[root->points_index_[0]].x,
                              Point.y - (*db)[root->points_index_[0]].y,
                              Point.z - (*db)[root->points_index_[0]].z);
       float dis = radius.squaredNorm();
       a.add_point(dis,root->points_index_[0]);
       //a.worst_dis = a.worst_dis < dis? a.worst_dis:dis;
       //判断现在的球是否在立方体内，如果在可以提前终止
       bool q = inside(root,Point,a.worst_dis);
      // cout << a.worst_dis_cap[0].dist << endl;
       return q;
    }
    //判断目标点所属象限
    int Coordinate = 0;
    if(Point.x > root->center_.x){
        Coordinate = Coordinate | 1;
    }
    if(Point.y > root->center_.y){
        Coordinate = Coordinate | 2;
    }
    if(Point.z > root->center_.z){
        Coordinate = Coordinate | 4;
    }
    //迭代寻找新的子象限
    if(octree_knn_search(root->child_root_[Coordinate],db,Point,a)) return true;
    //当发现最近的子象限都不能完全包裹最坏距离，那么就要扫描其他的子象限
    for(int i = 0;i<8;i++){
        //先排除刚才已经扫描过的象限
        if(i == Coordinate || root->child_root_[i] == nullptr) continue;
        //再排除球与立方体不相交的情况
        //cout << i << endl;
        if(false == overlap(root->child_root_[i],Point,a.worst_dis)) continue;
        //最后对这个象限进行计算worst_dis
        if(octree_knn_search(root->child_root_[i],db,Point,a)) return true;
    }

    //再次判断现在的球是否在立方体内，如果在可以提前终止
    return inside(root,Point,a.worst_dis);
}
int main(int argc, char const *argv[])
{
    
    ifstream fin;
    fin.open(db_list,ios::binary);
    if(!fin){
		cout<<"open error!"<<endl;
		return -1;
	}
	pcl::PointCloud<PointXYZ>::Ptr points (new pcl::PointCloud<PointXYZ>);
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_filtered(new pcl::PointCloud<pcl::PointXYZ>);
    pcl::VoxelGrid<pcl::PointXYZ> sor;
	int i;
    float x_min,x_max,y_min,y_max,z_min,z_max;
    x_min = com_float;
    x_max = -com_float;
    y_min = com_float;
    y_max = -com_float;
    z_min = com_float;
    z_max = -com_float;
	for (i=0; fin.good() && !fin.eof(); i++) {
		PointXYZ point;
		fin.read((char *) &point.x, 3*sizeof(float));
		points->push_back(point);
        x_min = (x_min < point.x)?x_min:point.x;
        x_max = (x_max > point.x)?x_max:point.x;
        y_min = (y_min < point.y)?y_min:point.y;
        y_max = (y_max > point.y)?y_max:point.y;
        z_min = (z_min < point.z)?z_min:point.z;
        z_max = (z_max > point.z)?z_max:point.z;
	}
        cout << "点云x范围=【" << x_min << ","
             << " "<< x_max << "】"
             << "点云y范围=【" << y_min << ","
             << " " << y_max << "】"
             << "点云z范围=【" << z_min << ","
             << " " << z_max << "】"
             << endl;
    /*降采样*/
    sor.setInputCloud(points);
    sor.setLeafSize(0.7f, 0.7f, 0.7f);
    sor.filter(*cloud_filtered);
    std::vector<int> points_index(cloud_filtered->size());
    points_index[0] = 0;
    std::partial_sum(points_index.begin(), points_index.end(), points_index.begin(), [](const int&a, int b) {return a + 1;});
    /*计算包含所有点的大立方体*/
    PointXYZ center((x_min+x_max)/2.,(y_min+y_max)/2.,(z_min+z_max)/2.);
    float extent = 0;
    extent = (y_max-y_min)<(x_max-x_min)?(x_max-x_min):(y_max-y_min);
    extent = extent<(z_max-z_min)?(z_max-z_min):extent;
    extent = ceil(extent / 2.);
    Octant* root;
    PointXYZ goals_point(21,0.1,0.63);//目标点
    int indexnumber = 0;
    cout << "输入需要的搜寻最近邻点的数量=" ;
    cin >> indexnumber;
    /*暴力搜索knn*/
    map<float,int> distance;
    chrono::steady_clock::time_point t5 = chrono::steady_clock::now();
    for(int k:points_index){
        Eigen::Vector3d radius(goals_point.x - (*cloud_filtered)[k].x,
                               goals_point.y - (*cloud_filtered)[k].y,
                               goals_point.z - (*cloud_filtered)[k].z);
        distance.insert(make_pair(radius.squaredNorm(),k));
    
    }
    cout << "最近的点的index:";
    int num=0;
  	for(map<float,int>::iterator it=distance.begin();it!=distance.end();++it){
		if(num++ == indexnumber) break;
        cout<<it->second<<" ";
    }
    cout << endl;
    chrono::steady_clock::time_point t6 = chrono::steady_clock::now();
    chrono::duration<double> time_used2 = chrono::duration_cast<chrono::duration<double>>(t6 - t5)*1000;
    cout << "暴力搜索用时 = " << time_used2.count() << " ms.    " <<endl;
    /*建立八叉树*/
    chrono::steady_clock::time_point t1 = chrono::steady_clock::now();
    root = build_octree(root,cloud_filtered,center,extent,points_index,0,-1);
    chrono::steady_clock::time_point t2 = chrono::steady_clock::now();
    chrono::duration<double> time_used = chrono::duration_cast<chrono::duration<double>>(t2 - t1)*1000;
    cout << "建立八叉树用时 = " << time_used.count() << " ms.    " << endl;
    result a(2e6,indexnumber);
    chrono::steady_clock::time_point t3 = chrono::steady_clock::now();
    octree_knn_search(root,cloud_filtered,goals_point,a);
    cout << "最近的点的index:";
    for(int i=0;i<a.size;i++){
        cout << a.worst_dis_cap[i].index  << " ";
    }
    cout << endl;
    chrono::steady_clock::time_point t4 = chrono::steady_clock::now();
    chrono::duration<double> time_used1 = chrono::duration_cast<chrono::duration<double>>(t4 - t3)*1000;
    cout << "knn用时 = " << time_used1.count() << " ms. " << endl;
    /*显示*/
    /* pcl::visualization::PCLVisualizer viewer("demo");
    int v1(0);
    viewer.createViewPort(0.0, 0.0, 0.5, 1.0, v1);
    pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> cloud_in_color_h(cloud_filtered, (int)0, (int)255 , (int)255 );
    viewer.addPointCloud(cloud_filtered, cloud_in_color_h, "cloud_in_v1", v1); 
    viewer.setSize(1280, 1024);

    while (!viewer.wasStopped())
    {
        cout << points->size() << " " << cloud_filtered->size()  << endl;
        viewer.spinOnce();
    }*/

    return 0;
}
