#include "io/config_xml.h"
#include "shape/variant.h"
#include "particle.h"
#include "configuration.h"
#include "BoundaryCondition.h"
#include "obj_loader.h"
#include "tinyxml2/tinyxml2.h"
#include <unistd.h>

class ShapePrintVisitor: public boost::static_visitor<std::string>{
public:
    std::string operator()(const shape::Polyhedron& poly)const{
        const char* source = poly.get_source();
        std::string ret;
        if(source){
            ret = "<shape type=\"polyhedron\" src=\"" + std::string(source) + "\"/>\n";
        }
        else{
            ret = "<shape type=\"polyhedron\">\n";
            ret += "<vertices>\n";
            for(auto vert: poly.vertices()){
                ret += "<vertex><x>" + std::to_string(vert[0]) + "</x><y>" + std::to_string(vert[1]) + "</y><z>" + std::to_string(vert[2]) + "</z></vertex>\n";
            }
            ret += "</vertices>\n";
            ret += "<faces>\n";
            for(auto face: poly.faces()){
                ret += "<face>\n";
                for(auto face_idx: face){
                    ret += "<fi>" + std::to_string(face_idx) + "</fi>\n";
                }
                ret += "</face>\n";
            }
            ret += "</faces>\n</shape>\n";
        }
        return ret;
    }

    std::string operator()(const shape::Sphere& sph)const{
        std::string ret("<shape type=\"sphere\">\n");
        ret += "<radius>" + std::to_string(sph.radius()) + "</radius>";
        ret += "</shape>\n";
        return ret;
    }

    std::string operator()(const shape::Box& box)const{
        std::string ret("<shape type=\"box\">\n");
        ret += "<x>" + std::to_string(2.0 * box.extent()[0]) + "</x>";
        ret += "<y>" + std::to_string(2.0 * box.extent()[1]) + "</y>";
        ret += "<z>" + std::to_string(2.0 * box.extent()[2]) + "</z>\n";
        ret += "</shape>\n";
        return ret;
    }

    std::string operator()(const shape::Cone& cone)const{
        std::string ret("<shape type=\"cone\">\n");
        ret += "<base_radius>" + std::to_string(cone.base_radius()) + "</base_radius>\n";
        ret += "<height>" + std::to_string(cone.height()) + "</height>\n";
        ret += "</shape>\n";
        return ret;
    }

    std::string operator()(const shape::Cylinder& cyl)const{
        std::string ret("<shape type=\"cylinder\">\n");
        ret += "<base_radius>" + std::to_string(cyl.base_radius()) + "</base_radius>\n";
        ret += "<height>" + std::to_string(cyl.height()) + "</height>\n";
        ret += "</shape>\n";
        return ret;
    }
};

void print_shape(FILE* fp, const shape::Variant& shape){
    auto str = boost::apply_visitor(ShapePrintVisitor(), shape);
    fprintf(fp, str.c_str());
}

void print_particle(FILE* fp, const Particle& particle, bool should_print_shape_id){
    fprintf(fp, "<particle>\n");

    fprintf(fp, "<pos>\n");
    fprintf(fp, "<x>%lf</x><y>%lf</y><z>%lf</z>\n", particle.xform.pos_[0], particle.xform.pos_[1], particle.xform.pos_[2]);
    fprintf(fp, "</pos>\n");

    //TODO: Perhaps skip this if the particle shape is a sphere?
    fprintf(fp, "<quat>\n");
    fprintf(fp, "<x>%lf</x><y>%lf</y><z>%lf</z><w>%lf</w>\n", particle.xform.rot_[0], particle.xform.rot_[1], particle.xform.rot_[2], particle.xform.rot_[3]);
    fprintf(fp, "</quat>\n");

    if(particle.xform.size_ != 1.0){
        fprintf(fp, "<size>%f</size>\n", particle.xform.size_);
    }

    if(should_print_shape_id) fprintf(fp, "<shape id=\"%lu\"/>\n", particle.shape_id);

    fprintf(fp, "</particle>\n");
}

bool xml_save_config(const char* filename, const Configuration& config){
    FILE* fp = fopen(filename, "w");
    if(!fp) return false;

    fprintf(fp, "<sim id=\"%d\">\n", getpid());{
        fprintf(fp, "<definitions>\n");{
            fprintf(fp, "<box shape=\"%s\" pbc=\"%s\">\n", "Rectangular", (false)? "false": "true");{
                const clam::Vec3d& size = config.pbc_.getSize();
                fprintf(fp, "<size>\n");
                fprintf(fp, "<x>%lf</x>\n", size[0]);
                fprintf(fp, "<y>%lf</y>\n", size[1]);
                fprintf(fp, "<z>%lf</z>\n", size[2]);
                fprintf(fp, "</size>\n");
            }fprintf(fp, "</box>\n");
            fprintf(fp, "<shapes>\n");{
                for(auto shape: config.shapes_){
                    print_shape(fp, *shape);
                }
            }fprintf(fp, "</shapes>\n");
        }fprintf(fp, "</definitions>\n");
        fprintf(fp, "<particles>\n");{
            for(auto particle: config.particles_){
                print_particle(fp, particle, config.shapes_.size() > 1);
            }
        }fprintf(fp, "</particles>\n");
    }fprintf(fp, "</sim>");

    fclose(fp);

    return true;
}

bool xml_load_config(const char* filename, Configuration& config){
    using namespace tinyxml2;
    XMLDocument doc;
    auto error = doc.LoadFile(filename);
    if(error != XML_NO_ERROR){
        printf("Error loading xml.\n");
        return false;
    }

    XMLElement* root = doc.RootElement();
    if(!root) return false;

    const XMLElement* definitions = root->FirstChildElement("definitions");
    if(!definitions) return false;
    //Load box
    {
        const XMLElement* box = definitions->FirstChildElement("box");
        if(!box) return false;

        //bool box_has_pbc = box->BoolAttribute("pbc");

        if(box->Attribute("shape", "Rectangular")){
            clam::Vec3d box_size(1.0);
            const XMLElement* size_node = box->FirstChildElement("size");
            size_node->FirstChildElement("x")->QueryDoubleText(&box_size[0]);
            size_node->FirstChildElement("y")->QueryDoubleText(&box_size[1]);
            size_node->FirstChildElement("z")->QueryDoubleText(&box_size[2]);
            config.pbc_.setSize(box_size);
        }
        else{
            printf("Unknown box type: %s.\n", box->Attribute("shape"));
            return false;
        }
    }
    //Load shapes
    {
        const XMLElement* shapes = definitions->FirstChildElement("shapes");
        if(!shapes) return false;

        for(const XMLElement* shape = shapes->FirstChildElement(); shape != nullptr; shape = shape->NextSiblingElement()){
            if(shape->Attribute("type", "polyhedron")){
                std::vector<clam::Vec3d> vertices;
                std::vector<std::vector<unsigned int>> faces;

                if(const char* source = shape->Attribute("src")){
                    if(!load_obj(source, vertices, faces)) return false;
                }
                else{
                    const XMLElement* v_node = shape->FirstChildElement("vertices");
                    if(!v_node) return false;

                    for(const XMLElement* vert = v_node->FirstChildElement();
                        vert != nullptr;
                        vert = vert->NextSiblingElement()
                    ){
                        clam::Vec3d vertex;
                        vert->FirstChildElement("x")->QueryDoubleText(&vertex[0]);
                        vert->FirstChildElement("y")->QueryDoubleText(&vertex[1]);
                        vert->FirstChildElement("z")->QueryDoubleText(&vertex[2]);
                        vertices.push_back(vertex);
                    }

                    const XMLElement* f_node = shape->FirstChildElement("faces");
                    if(!f_node) return false;

                    //Loop over faces
                    for(const XMLElement* fi_node = f_node->FirstChildElement();
                        fi_node != nullptr;
                        fi_node = fi_node->NextSiblingElement()
                    ){
                        std::vector<unsigned int> face;
                        //Loop over fi
                        for(const XMLElement* fi = fi_node->FirstChildElement(); fi != nullptr; fi = fi->NextSiblingElement()){
                            unsigned int f_idx = 0;
                            fi->QueryUnsignedText(&f_idx);
                            face.push_back(f_idx);
                        }
                        faces.push_back(face);
                    }
                }

                config.shapes_.push_back(new shape::Variant(shape::Polyhedron(vertices, faces)));
            }
            else if(shape->Attribute("type", "sphere")){
                config.shapes_.push_back(new shape::Variant(shape::Sphere()));
            }
            else if(shape->Attribute("type", "box")){
                clam::Vec3d size(0.0);
                shape->FirstChildElement("x")->QueryDoubleText(&size[0]);
                shape->FirstChildElement("y")->QueryDoubleText(&size[1]);
                shape->FirstChildElement("z")->QueryDoubleText(&size[2]);
                config.shapes_.push_back(new shape::Variant(shape::Box(size)));
            }
            else if(shape->Attribute("type", "cone")){
                double base_radius, height;
                shape->FirstChildElement("base_radius")->QueryDoubleText(&base_radius);
                shape->FirstChildElement("height")->QueryDoubleText(&height);
                config.shapes_.push_back(new shape::Variant(shape::Cone(base_radius, height)));
            }
            else if(shape->Attribute("type", "cylinder")){
                double base_radius, height;
                shape->FirstChildElement("base_radius")->QueryDoubleText(&base_radius);
                shape->FirstChildElement("height")->QueryDoubleText(&height);
                config.shapes_.push_back(new shape::Variant(shape::Cylinder(base_radius, height)));
            }
            else{
                printf("Unknown shape type: %s.\n", shape->Attribute("type"));
                return false;
            }
        }
    }

    //Load particles
    {
        const XMLElement* particles = root->FirstChildElement("particles");
        if(!particles) return false;

        int i = 0;
        for(const XMLElement* particle_element = particles->FirstChildElement();
            particle_element != nullptr;
            particle_element = particle_element->NextSiblingElement(), ++i
        ){
            Particle p;
            double x, y, z, w;

            const XMLElement* pos_element = particle_element->FirstChildElement("pos");
            pos_element->FirstChildElement("x")->QueryDoubleText(&x);
            pos_element->FirstChildElement("y")->QueryDoubleText(&y);
            pos_element->FirstChildElement("z")->QueryDoubleText(&z);
            p.xform.pos_ = clam::Vec3d(x, y, z);

            const XMLElement* rot_element = particle_element->FirstChildElement("quat");
            rot_element->FirstChildElement("x")->QueryDoubleText(&x);
            rot_element->FirstChildElement("y")->QueryDoubleText(&y);
            rot_element->FirstChildElement("z")->QueryDoubleText(&z);
            rot_element->FirstChildElement("w")->QueryDoubleText(&w);
            p.xform.rot_ = clam::Quatd(x, y, z, w);

            if(config.shapes_.size() > 1){
                p.shape_id = particle_element->FirstChildElement("shape")->UnsignedAttribute("id");
            }
            else{
                p.shape_id = 0;
            }

            const XMLElement* size_element = particle_element->FirstChildElement("size");
            if(!size_element) p.xform.size_ = 1.0;
            else size_element->QueryDoubleText(&p.xform.size_);

            config.particles_.push_back(p);
        }
    }

    return true;
}

