#pragma once
#include "igl/opengl/glfw/Display.h"

static void glfw_mouse_press(GLFWwindow* window, int button, int action, int modifier)
{

  Renderer* rndr = (Renderer*) glfwGetWindowUserPointer(window);
  
  if (action == GLFW_PRESS)
  {
	  double x2, y2;
	  glfwGetCursorPos(window, &x2, &y2);
	  igl::opengl::glfw::Viewer* scn = rndr->GetScene();
	  bool found = false , first_object = true;
	  int i = 0, savedIndx = scn->selected_data_index, closest_index = scn->selected_data_index;
	  float new_depth,old_depth = 0;
	  for (; i < scn->data_list.size() ;i++)
	  {
		  scn->selected_data_index = i;
		  if (rndr->Picking(x2, y2)) {
			  found = true;
			  new_depth = rndr->get_currently_picked_object_z(x2, y2); //added function.
			  if (first_object) {
				  old_depth = new_depth;
				  closest_index = i;
				  first_object = false;
			  }
			  else if (old_depth < new_depth) {
				  closest_index = i;
				  old_depth = new_depth;
			  }
		  }
	  }
	  if (!found) {
		  std::cout << "scene selected" << std::endl;
		  scn->selected_data_index = savedIndx;
		  scn->picked = true;
	  }
	  else {
		  std::cout << "found " << (scn->selected_data_index = closest_index) << std::endl;
		  scn->picked = false;
	  }
	  rndr->UpdatePosition(x2, y2);
  }
}


 void glfw_mouse_move(GLFWwindow* window, double x, double y)
{
	 Renderer* rndr = (Renderer*)glfwGetWindowUserPointer(window);
	 rndr->UpdatePosition(x, y);
	 if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS)
	 {
		 rndr->MouseProcessing(GLFW_MOUSE_BUTTON_RIGHT);
	 }
	 else if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS)
	 {
		 rndr->MouseProcessing(GLFW_MOUSE_BUTTON_LEFT);
	 }
}

static void glfw_mouse_scroll(GLFWwindow* window, double x, double y)
{
	Renderer* rndr = (Renderer*)glfwGetWindowUserPointer(window);
	if (rndr->GetScene()->assignment3 || rndr->GetScene()->picked) rndr->GetScene()->MyTranslate(Eigen::Vector3f(0, 0, y * 0.1), true);
	else rndr->GetScene()->data().MyScale(Eigen::Vector3f(1 + y * 0.01,1 + y * 0.01,1+y*0.01));

}

void glfw_window_size(GLFWwindow* window, int width, int height)
{
	Renderer* rndr = (Renderer*)glfwGetWindowUserPointer(window);
  rndr->post_resize(window,width, height);

}


static void glfw_key_callback(GLFWwindow* window, int key, int scancode, int action, int modifier)
{
	Renderer* rndr = (Renderer*) glfwGetWindowUserPointer(window);
	igl::opengl::glfw::Viewer* scn = rndr->GetScene();
	if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
		glfwSetWindowShouldClose(window, GL_TRUE);

	else if(action == GLFW_PRESS || action == GLFW_REPEAT)
		switch (key)
		{
		case 'A':
		case 'a':
		{
			rndr->core().is_animating = !rndr->core().is_animating;
			break;
		}
		case 'F':
		case 'f':
		{
			scn->data().set_face_based(!scn->data().face_based);
			break;
		}
		case 'I':
		case 'i':
		{
			scn->data().dirty |= igl::opengl::MeshGL::DIRTY_NORMAL;
			scn->data().invert_normals = !scn->data().invert_normals;
			break;
		}
		case 'L':
		case 'l':
		{
			rndr->core().toggle(scn->data().show_lines);
			break;
		}
		case 'O':
		case 'o':
		{
			rndr->core().orthographic = !rndr->core().orthographic;
			break;
		}
		case '1':
		case '2':
		{
			scn->selected_data_index =
				(scn->selected_data_index + scn->data_list.size() + (key == '2' ? 1 : -1)) % scn->data_list.size();
			scn->picked = false;
			break;
		}
		case '[':
		case ']':
		{
			rndr->ChangeCamera(key);
			break;
		}
		case ';':
			scn->data().show_vertid = !scn->data().show_vertid;
			break;
		case ':':
			scn->data().show_faceid = !scn->data().show_faceid;
			break;
		case 's':
			if(scn->enable_semplefication) scn->simpleficate_mesh(std::ceil(0.05 * scn->data().Q->size()));
			else std::cout << "simplefication disabled" << std::endl;
			break;
		case ' ':
			scn->animating_ik = !scn->animating_ik;
			break;
		case 'p':
		case 'P':
			std::cout << "the rotation matrix of link " << scn->selected_data_index << " is :" << std::endl;
			std::cout << scn->get_parent_link_Rot(scn->selected_data_index) * scn->data().GetRotationTrans() << std::endl;
			break;
		case 't':
		case 'T':
			for (int link = 0; link < scn->link_number; link++) {
				std::cout
					<< "------------ Tip " << link << " location ------------------" << std::endl
					<< (scn->MakeTrans() * scn->get_parent_link_T(link) * scn->data_list[link].MakeTrans() * Eigen::Vector4f(0, 0.8, 0, 1)).block<3,1>(0,0).transpose()
					<< std::endl;
			}
			break;
		case 'q':
		case 'Q':
			scn->data_list[3].MyTranslate((scn->MakeTrans() * scn->data_list[1].MakeTrans()).block<3, 3>(0, 0) *Eigen::Vector3f(0, -0.8, 0),true);
			break;
		case 'd':
		case 'D':
			std::cout
				<< "------------ Destination location ------------------" << std::endl
				<< (scn->data_list[scn->sphere_index].MakeTrans() * Eigen::Vector4f(0, 0, 0, 1)).block<3, 1>(0, 0).transpose()
				<< std::endl;
		break;
		case 263:
			scn->data().MyRotate(Eigen::Vector3f(0, 1, 0), -0.1);
			break;
		case 262:
			scn->data().MyRotate(Eigen::Vector3f(0, 1, 0), 0.1);
			break;
		case 265:
			scn->data().MyRotate(Eigen::Vector3f(1, 0, 0), -0.1);
			break;
		case 264:
			scn->data().MyRotate(Eigen::Vector3f(1, 0, 0),0.1);
			break;
		default:break;//do nothing
		}
}


void Init(Display& display)
{
	display.AddKeyCallBack(glfw_key_callback);
	display.AddMouseCallBacks(glfw_mouse_press, glfw_mouse_scroll, glfw_mouse_move);
	display.AddResizeCallBack(glfw_window_size);
}
