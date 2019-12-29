
#include "igl/opengl/glfw/renderer.h"
#include "tutorial/sandBox/inputManager.h"


int main(int argc, char *argv[])
{
  Display *disp = new Display(1000, 800, "Wellcome");
  Renderer renderer;
  igl::opengl::glfw::Viewer viewer;
  //-----------------------------------
  bool enable_simplefication = false; // enables and disables simplefication option
  bool render_arm_ass3 = true;
  viewer.load_mesh_from_configuration(enable_simplefication , render_arm_ass3);
  //-----------------------------------
  Init(*disp);
  renderer.init(&viewer);
  disp->SetRenderer(&renderer);
  disp->launch_rendering(true);
  
  delete disp;
}
