/************************************************************************************

Authors     :   Bradley Austin Davis <bdavis@saintandreas.org>
Copyright   :   Copyright Brad Davis. All Rights reserved.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

************************************************************************************/

#include <iostream>
#include <memory>
#include <exception>
#include <algorithm>

#include <Windows.h>

#define __STDC_FORMAT_MACROS 1

#define FAIL(X) throw std::runtime_error(X)

///////////////////////////////////////////////////////////////////////////////
//
// GLM is a C++ math library meant to mirror the syntax of GLSL 
//

#include <glm/glm.hpp>
#include <glm/gtc/noise.hpp>
#include <glm/gtx/rotate_vector.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/quaternion.hpp>
#include <glm/gtx/quaternion.hpp>
#include "Skybox.h"

// Import the most commonly used types into the default namespace
using glm::ivec3;
using glm::ivec2;
using glm::uvec2;
using glm::mat3;
using glm::mat4;
using glm::vec2;
using glm::vec3;
using glm::vec4;
using glm::quat;

///////////////////////////////////////////////////////////////////////////////
//
// GLEW gives cross platform access to OpenGL 3.x+ functionality.  
//

#include <GL/glew.h>

bool checkFramebufferStatus(GLenum target = GL_FRAMEBUFFER)
{
  GLuint status = glCheckFramebufferStatus(target);
  switch (status)
  {
  case GL_FRAMEBUFFER_COMPLETE:
    return true;
    break;

  case GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT:
    std::cerr << "framebuffer incomplete attachment" << std::endl;
    break;

  case GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT:
    std::cerr << "framebuffer missing attachment" << std::endl;
    break;

  case GL_FRAMEBUFFER_INCOMPLETE_DRAW_BUFFER:
    std::cerr << "framebuffer incomplete draw buffer" << std::endl;
    break;

  case GL_FRAMEBUFFER_INCOMPLETE_READ_BUFFER:
    std::cerr << "framebuffer incomplete read buffer" << std::endl;
    break;

  case GL_FRAMEBUFFER_INCOMPLETE_MULTISAMPLE:
    std::cerr << "framebuffer incomplete multisample" << std::endl;
    break;

  case GL_FRAMEBUFFER_INCOMPLETE_LAYER_TARGETS:
    std::cerr << "framebuffer incomplete layer targets" << std::endl;
    break;

  case GL_FRAMEBUFFER_UNSUPPORTED:
    std::cerr << "framebuffer unsupported internal format or image" << std::endl;
    break;

  default:
    std::cerr << "other framebuffer error" << std::endl;
    break;
  }

  return false;
}

bool checkGlError()
{
  GLenum error = glGetError();
  if (!error)
  {
    return false;
  }
  else
  {
    switch (error)
    {
    case GL_INVALID_ENUM:
      std::cerr <<
        ": An unacceptable value is specified for an enumerated argument.The offending command is ignored and has no other side effect than to set the error flag.";
      break;
    case GL_INVALID_VALUE:
      std::cerr <<
        ": A numeric argument is out of range.The offending command is ignored and has no other side effect than to set the error flag";
      break;
    case GL_INVALID_OPERATION:
      std::cerr <<
        ": The specified operation is not allowed in the current state.The offending command is ignored and has no other side effect than to set the error flag..";
      break;
    case GL_INVALID_FRAMEBUFFER_OPERATION:
      std::cerr <<
        ": The framebuffer object is not complete.The offending command is ignored and has no other side effect than to set the error flag.";
      break;
    case GL_OUT_OF_MEMORY:
      std::cerr <<
        ": There is not enough memory left to execute the command.The state of the GL is undefined, except for the state of the error flags, after this error is recorded.";
      break;
    case GL_STACK_UNDERFLOW:
      std::cerr <<
        ": An attempt has been made to perform an operation that would cause an internal stack to underflow.";
      break;
    case GL_STACK_OVERFLOW:
      std::cerr << ": An attempt has been made to perform an operation that would cause an internal stack to overflow.";
      break;
    }
    return true;
  }
}

void glDebugCallbackHandler(GLenum source, GLenum type, GLuint id, GLenum severity, GLsizei length, const GLchar* msg,
                            GLvoid* data)
{
  OutputDebugStringA(msg);
  std::cout << "debug call: " << msg << std::endl;
}

//////////////////////////////////////////////////////////////////////
//
// GLFW provides cross platform window creation
//

#include <GLFW/glfw3.h>

namespace glfw
{
  inline GLFWwindow* createWindow(const uvec2& size, const ivec2& position = ivec2(INT_MIN))
  {
    GLFWwindow* window = glfwCreateWindow(size.x, size.y, "glfw", nullptr, nullptr);
    if (!window)
    {
      FAIL("Unable to create rendering window");
    }
    if ((position.x > INT_MIN) && (position.y > INT_MIN))
    {
      glfwSetWindowPos(window, position.x, position.y);
    }
    return window;
  }
}

// A class to encapsulate using GLFW to handle input and render a scene
class GlfwApp
{
protected:
  uvec2 windowSize;
  ivec2 windowPosition;
  GLFWwindow* window{nullptr};
  unsigned int frame{0};

public:
  GlfwApp()
  {
    // Initialize the GLFW system for creating and positioning windows
    if (!glfwInit())
    {
      FAIL("Failed to initialize GLFW");
    }
    glfwSetErrorCallback(ErrorCallback);
  }

  virtual ~GlfwApp()
  {
    if (nullptr != window)
    {
      glfwDestroyWindow(window);
    }
    glfwTerminate();
  }

  virtual int run()
  {
    preCreate();

    window = createRenderingTarget(windowSize, windowPosition);

    if (!window)
    {
      std::cout << "Unable to create OpenGL window" << std::endl;
      return -1;
    }

    postCreate();

    initGl();

    while (!glfwWindowShouldClose(window))
    {
      ++frame;
      glfwPollEvents();
      update();
      draw();
      finishFrame();
    }

    shutdownGl();

    return 0;
  }

protected:
  virtual GLFWwindow* createRenderingTarget(uvec2& size, ivec2& pos) = 0;

  virtual void draw() = 0;

  void preCreate()
  {
    glfwWindowHint(GLFW_DEPTH_BITS, 16);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_DEBUG_CONTEXT, true);
  }

  void postCreate()
  {
    glfwSetWindowUserPointer(window, this);
    glfwSetKeyCallback(window, KeyCallback);
    glfwSetMouseButtonCallback(window, MouseButtonCallback);
    glfwMakeContextCurrent(window);

    // Initialize the OpenGL bindings
    // For some reason we have to set this experminetal flag to properly
    // init GLEW if we use a core context.
    glewExperimental = GL_TRUE;
    if (0 != glewInit())
    {
      FAIL("Failed to initialize GLEW");
    }
    glGetError();

    if (GLEW_KHR_debug)
    {
      GLint v;
      glGetIntegerv(GL_CONTEXT_FLAGS, &v);
      if (v & GL_CONTEXT_FLAG_DEBUG_BIT)
      {
        //glDebugMessageCallback(glDebugCallbackHandler, this);
      }
    }
  }

  virtual void initGl()
  {
  }

  virtual void shutdownGl()
  {
  }

  virtual void finishFrame()
  {
    glfwSwapBuffers(window);
  }

  virtual void destroyWindow()
  {
    glfwSetKeyCallback(window, nullptr);
    glfwSetMouseButtonCallback(window, nullptr);
    glfwDestroyWindow(window);
  }

  virtual void onKey(int key, int scancode, int action, int mods)
  {
    if (GLFW_PRESS != action)
    {
      return;
    }

    switch (key)
    {
    case GLFW_KEY_ESCAPE:
      glfwSetWindowShouldClose(window, 1);
      return;
    }
  }

  virtual void update()
  {
  }

  virtual void onMouseButton(int button, int action, int mods)
  {
  }

protected:
  virtual void viewport(const ivec2& pos, const uvec2& size)
  {
    glViewport(pos.x, pos.y, size.x, size.y);
  }

private:

  static void KeyCallback(GLFWwindow* window, int key, int scancode, int action, int mods)
  {
    GlfwApp* instance = (GlfwApp *)glfwGetWindowUserPointer(window);
    instance->onKey(key, scancode, action, mods);
  }

  static void MouseButtonCallback(GLFWwindow* window, int button, int action, int mods)
  {
    GlfwApp* instance = (GlfwApp *)glfwGetWindowUserPointer(window);
    instance->onMouseButton(button, action, mods);
  }

  static void ErrorCallback(int error, const char* description)
  {
    FAIL(description);
  }
};

//////////////////////////////////////////////////////////////////////
//
// The Oculus VR C API provides access to information about the HMD
//


#include <OVR_CAPI.h>
#include <OVR_CAPI_GL.h>

// Head Tracking Modes (BUTTON B)
enum HeadTrackingMode { REGULAR, ORIENTATION_ONLY, POSITION_ONLY, NO_TRACKING };
HeadTrackingMode headMode;

namespace ovr
{
  // Convenience method for looping over each eye with a lambda
  template <typename Function>
  inline void for_each_eye(Function function)
  {
    for (ovrEyeType eye = ovrEyeType::ovrEye_Left;
         eye < ovrEyeType::ovrEye_Count;
         eye = static_cast<ovrEyeType>(eye + 1))
    {
      function(eye);
    }
  }

  inline mat4 toGlm(const ovrMatrix4f& om)
  {
    return glm::transpose(glm::make_mat4(&om.M[0][0]));
  }

  inline mat4 toGlm(const ovrFovPort& fovport, float nearPlane = 0.01f, float farPlane = 10000.0f)
  {
    return toGlm(ovrMatrix4f_Projection(fovport, nearPlane, farPlane, true));
  }

  inline vec3 toGlm(const ovrVector3f& ov)
  {
    return glm::make_vec3(&ov.x);
  }

  inline vec2 toGlm(const ovrVector2f& ov)
  {
    return glm::make_vec2(&ov.x);
  }

  inline uvec2 toGlm(const ovrSizei& ov)
  {
    return uvec2(ov.w, ov.h);
  }

  inline quat toGlm(const ovrQuatf& oq)
  {
    return glm::make_quat(&oq.x);
  }

  inline mat4 toGlm(const ovrPosef& op)
  {
    mat4 orientation = glm::mat4_cast(toGlm(op.Orientation));
    mat4 translation = glm::translate(mat4(), ovr::toGlm(op.Position));
	
	switch (headMode) {
	case REGULAR:
		break;
	case ORIENTATION_ONLY:
		translation = glm::mat4(1.0f);
		break;
	case POSITION_ONLY:
		orientation = glm::mat4(1.0f);
		break;
	case NO_TRACKING:
		translation = glm::mat4(1.0f);
		orientation = glm::mat4(1.0f);
		break;
	}
    return translation * orientation;
	
  }

  inline ovrMatrix4f fromGlm(const mat4& m)
  {
    ovrMatrix4f result;
    mat4 transposed(glm::transpose(m));
    memcpy(result.M, &(transposed[0][0]), sizeof(float) * 16);
    return result;
  }

  inline ovrVector3f fromGlm(const vec3& v)
  {
    ovrVector3f result;
    result.x = v.x;
    result.y = v.y;
    result.z = v.z;
    return result;
  }

  inline ovrVector2f fromGlm(const vec2& v)
  {
    ovrVector2f result;
    result.x = v.x;
    result.y = v.y;
    return result;
  }

  inline ovrSizei fromGlm(const uvec2& v)
  {
    ovrSizei result;
    result.w = v.x;
    result.h = v.y;
    return result;
  }

  inline ovrQuatf fromGlm(const quat& q)
  {
    ovrQuatf result;
    result.x = q.x;
    result.y = q.y;
    result.z = q.z;
    result.w = q.w;
    return result;
  }
}

#include <queue>

class RiftManagerApp
{
protected:
  ovrSession _session;
  ovrHmdDesc _hmdDesc;
  ovrGraphicsLuid _luid;

public:
  RiftManagerApp()
  {
    if (!OVR_SUCCESS(ovr_Create(&_session, &_luid)))
    {
      FAIL("Unable to create HMD session");
    }

    _hmdDesc = ovr_GetHmdDesc(_session);
  }

  ~RiftManagerApp()
  {
    ovr_Destroy(_session);
    _session = nullptr;
  }
};

// Hand Tracking
double displayMidpointSeconds;
ovrTrackingState trackState;
ovrInputState  inputState;

unsigned int handStatus[2]; //Hand Status
ovrPosef handPoses[2]; // Hand Poses
ovrVector3f handPosition[2];
ovrQuatf handRotation[2];


// Game Mode
enum GameMode { PREPARE, ON, END };
GameMode gameMode;

// Prepare Modes (BUTTON B, BUTTON A)
// Button B
enum DirectionMode { NONE, UP, RIGHT, DOWN, LEFT };
DirectionMode directionMode;
// Button A
enum WarshipMode { A_Ship, B_Ship, C_Ship, S_Ship, P_Ship };
WarshipMode warshipMode;
// BUTTON X
enum SelectingMode { SELECTING, DONE };
SelectingMode selectingMode;

// EYE Modes (BUTTON A)
enum SkyboxMode { SKYBOX_ENTIRE, SKYBOX_STEREO, SKYBOX_MONO, SKYBOX_MONO_CUSTOMIZED_1, SKYBOX_MONO_CUSTOMIZED_2 };
SkyboxMode skyboxMode;


// Scene Mode (BUTTON X)
enum ViewingMode { STEREO, MONO, LEFT_EYE_ONLY, RIGHT_EYE_ONLY, INVERTED_EYES };
ViewingMode viewMode;

class RiftApp : public GlfwApp, public RiftManagerApp
{
public:
	// IDO 
	float DAFAULT_IOD;
	const float MAX_IOD = 0.15f;
	const float MIN_IOD = -0.05f;
	float IOD;

	

	// Tracking Lag
	unsigned int Num_TrackingLagFrames;
	unsigned int NumTrackingFrames;
	std::queue<ovrPosef> LeftHandFrames;
	std::queue<ovrPosef> RightHandFrames;
	std::queue<ovrPosef> LeftEyeFrames;
	std::queue<ovrPosef> RightEyeFrames;

	// Rendering Delay
	unsigned int Num_RenderingDelayFrames;
	unsigned int NumRenderingFrames;

	// Smooth Mode
	unsigned int NumAverageMovingFrames;
	unsigned int NumSmoothFrames;
	std::vector<ovrPosef> LeftHandFrames_Smooth;
	std::vector<ovrPosef> RightHandFrames_Smooth;
	

	ovrPosef oldHandLeft, oldHandRight, currHandLeft, currHandRight;

private:
  GLuint _fbo{0};
  GLuint _depthBuffer{0};
  ovrTextureSwapChain _eyeTexture;

  GLuint _mirrorFbo{0};
  ovrMirrorTexture _mirrorTexture;

  ovrEyeRenderDesc _eyeRenderDescs[2];

  mat4 _eyeProjections[2];

  ovrLayerEyeFov _sceneLayer;
  ovrViewScaleDesc _viewScaleDesc;

  uvec2 _renderTargetSize;
  uvec2 _mirrorSize;

  ovrPosef oldEyeLeft, oldEyeRight, currEyeLeft, currEyeRight;
  
  

public:

  RiftApp()
  {
    using namespace ovr;
    _viewScaleDesc.HmdSpaceToWorldScaleInMeters = 1.0f;

    memset(&_sceneLayer, 0, sizeof(ovrLayerEyeFov));
    _sceneLayer.Header.Type = ovrLayerType_EyeFov;
    _sceneLayer.Header.Flags = ovrLayerFlag_TextureOriginAtBottomLeft;

    ovr::for_each_eye([&](ovrEyeType eye)
    {
      ovrEyeRenderDesc& erd = _eyeRenderDescs[eye] = ovr_GetRenderDesc(_session, eye, _hmdDesc.DefaultEyeFov[eye]);
      ovrMatrix4f ovrPerspectiveProjection =
        ovrMatrix4f_Projection(erd.Fov, 0.01f, 1000.0f, ovrProjection_ClipRangeOpenGL);
      _eyeProjections[eye] = ovr::toGlm(ovrPerspectiveProjection);
      _viewScaleDesc.HmdToEyePose[eye] = erd.HmdToEyePose;

      ovrFovPort& fov = _sceneLayer.Fov[eye] = _eyeRenderDescs[eye].Fov;
      auto eyeSize = ovr_GetFovTextureSize(_session, eye, fov, 1.0f);
      _sceneLayer.Viewport[eye].Size = eyeSize;
      _sceneLayer.Viewport[eye].Pos = {(int)_renderTargetSize.x, 0};

      _renderTargetSize.y = std::max(_renderTargetSize.y, (uint32_t)eyeSize.h);
      _renderTargetSize.x += eyeSize.w;
    });
    // Make the on screen window 1/4 the resolution of the render target
    _mirrorSize = _renderTargetSize;
    _mirrorSize /= 4;

	// IOD
	DAFAULT_IOD = std::abs(_viewScaleDesc.HmdToEyePose[0].Position.x - _viewScaleDesc.HmdToEyePose[1].Position.x);
	IOD = DAFAULT_IOD;

	// Tracking Lag
	Num_TrackingLagFrames = 0;
	NumTrackingFrames = 0;

	// Rendering Delay
	Num_RenderingDelayFrames = 0;
	NumRenderingFrames = 0;

	// Smooth Mode
	NumAverageMovingFrames = 1;
	NumSmoothFrames = 0;
  }

protected:
  GLFWwindow* createRenderingTarget(uvec2& outSize, ivec2& outPosition) override
  {
    return glfw::createWindow(_mirrorSize);
  }

  void initGl() override
  {
    GlfwApp::initGl();

    // Disable the v-sync for buffer swap
    glfwSwapInterval(0);

    ovrTextureSwapChainDesc desc = {};
    desc.Type = ovrTexture_2D;
    desc.ArraySize = 1;
    desc.Width = _renderTargetSize.x;
    desc.Height = _renderTargetSize.y;
    desc.MipLevels = 1;
    desc.Format = OVR_FORMAT_R8G8B8A8_UNORM_SRGB;
    desc.SampleCount = 1;
    desc.StaticImage = ovrFalse;
    ovrResult result = ovr_CreateTextureSwapChainGL(_session, &desc, &_eyeTexture);
    _sceneLayer.ColorTexture[0] = _eyeTexture;
    if (!OVR_SUCCESS(result))
    {
      FAIL("Failed to create swap textures");
    }

    int length = 0;
    result = ovr_GetTextureSwapChainLength(_session, _eyeTexture, &length);
    if (!OVR_SUCCESS(result) || !length)
    {
      FAIL("Unable to count swap chain textures");
    }
    for (int i = 0; i < length; ++i)
    {
      GLuint chainTexId;
      ovr_GetTextureSwapChainBufferGL(_session, _eyeTexture, i, &chainTexId);
      glBindTexture(GL_TEXTURE_2D, chainTexId);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    }
    glBindTexture(GL_TEXTURE_2D, 0);

    // Set up the framebuffer object
    glGenFramebuffers(1, &_fbo);
    glGenRenderbuffers(1, &_depthBuffer);
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, _fbo);
    glBindRenderbuffer(GL_RENDERBUFFER, _depthBuffer);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT16, _renderTargetSize.x, _renderTargetSize.y);
    glBindRenderbuffer(GL_RENDERBUFFER, 0);
    glFramebufferRenderbuffer(GL_DRAW_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, _depthBuffer);
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);

    ovrMirrorTextureDesc mirrorDesc;
    memset(&mirrorDesc, 0, sizeof(mirrorDesc));
    mirrorDesc.Format = OVR_FORMAT_R8G8B8A8_UNORM_SRGB;
    mirrorDesc.Width = _mirrorSize.x;
    mirrorDesc.Height = _mirrorSize.y;
    if (!OVR_SUCCESS(ovr_CreateMirrorTextureGL(_session, &mirrorDesc, &_mirrorTexture)))
    {
      FAIL("Could not create mirror texture");
    }
    glGenFramebuffers(1, &_mirrorFbo);
  }

  void onKey(int key, int scancode, int action, int mods) override
  {
    if (GLFW_PRESS == action)
      switch (key)
      {
      case GLFW_KEY_R:
        ovr_RecenterTrackingOrigin(_session);
        return;
      }

    GlfwApp::onKey(key, scancode, action, mods);
  }

  void draw() final override
  {
	 
	// HAND TRACKING
	displayMidpointSeconds = ovr_GetPredictedDisplayTime(_session, 0);
	trackState = ovr_GetTrackingState(_session, displayMidpointSeconds, ovrTrue);
	// Hand Status
	handStatus[0] = trackState.HandStatusFlags[0];
	handStatus[1] = trackState.HandStatusFlags[1];
	// Hand Poses
	handPoses[0] = trackState.HandPoses[0].ThePose;
	handPoses[1] = trackState.HandPoses[1].ThePose;
	// Hand Position
	handPosition[0] = handPoses[0].Position;
	handPosition[1] = handPoses[1].Position;
	// Hand Rotation
	handRotation[0] = handPoses[0].Orientation;
	handRotation[1] = handPoses[1].Orientation;

	// Eye Tracking
	ovrPosef eyePoses[2];
	ovrPosef originEyePoses[2];
	ovr_GetEyePoses(_session, frame, true, _viewScaleDesc.HmdToEyePose, eyePoses, &_sceneLayer.SensorSampleTime);
	originEyePoses[ovrEye_Left] = eyePoses[ovrEye_Left];
	originEyePoses[ovrEye_Right] = eyePoses[ovrEye_Right];

	// EXTRA CREDIT 2: Smooth Mode
	if (NumAverageMovingFrames > 1) {
		if (NumSmoothFrames == NumAverageMovingFrames) {
			NumSmoothFrames = 0;
			

			ovrVector3f LPos;
			ovrVector3f RPos;

			LPos = LeftHandFrames_Smooth[0].Position;
			RPos = RightHandFrames_Smooth[0].Position;
			for (unsigned int i = 1; i < LeftHandFrames_Smooth.size(); i++) {
				LPos.x += LeftHandFrames_Smooth[i].Position.x;
				LPos.y += LeftHandFrames_Smooth[i].Position.y;
				LPos.z += LeftHandFrames_Smooth[i].Position.z;
				RPos.x += RightHandFrames_Smooth[i].Position.x;
				RPos.y += RightHandFrames_Smooth[i].Position.y;
				RPos.z += RightHandFrames_Smooth[i].Position.z;
			}

			LPos.x = LPos.x / NumAverageMovingFrames;
			LPos.y = LPos.y / NumAverageMovingFrames;
			LPos.z = LPos.z / NumAverageMovingFrames;
			RPos.x = RPos.x / NumAverageMovingFrames;
			RPos.y = RPos.y / NumAverageMovingFrames;
			RPos.z = RPos.z / NumAverageMovingFrames;

			currHandLeft = { handPoses[ovrHand_Left].Orientation,LPos };
			currHandRight = { handPoses[ovrHand_Right].Orientation,RPos };
			oldHandLeft = { currHandLeft.Orientation, currHandLeft.Position };
			oldHandRight = { currHandRight.Orientation, currHandRight.Position };
			handPoses[ovrHand_Left].Position = LPos;
			handPoses[ovrHand_Right].Position = RPos;

			
			LeftHandFrames_Smooth.clear();
			RightHandFrames_Smooth.clear();
			//std::cout << "SMOOTH MODE" << std::endl;
		}
		else {
			NumSmoothFrames++;
			LeftHandFrames_Smooth.push_back(handPoses[ovrHand_Left]);
			RightHandFrames_Smooth.push_back(handPoses[ovrHand_Right]);
			//std::cout << LeftHandFrames_Smooth.size() << std::endl;
		}
	}

	// Tracking Lag
	if (Num_TrackingLagFrames > 0) {
		if (LeftHandFrames.size() < Num_TrackingLagFrames) {
			LeftHandFrames.push(handPoses[ovrHand_Left]);
			RightHandFrames.push(handPoses[ovrHand_Right]);
			LeftEyeFrames.push(eyePoses[ovrEye_Left]);
			RightEyeFrames.push(eyePoses[ovrEye_Right]);
		}
		else {
			currHandLeft = LeftHandFrames.front();
			LeftHandFrames.pop();
			LeftHandFrames.push(handPoses[ovrHand_Left]);
			handPoses[ovrHand_Left] = currHandLeft;

			currHandRight = RightHandFrames.front();
			RightHandFrames.pop();
			RightHandFrames.push(handPoses[ovrHand_Right]);
			handPoses[ovrHand_Right] = currHandRight;

			currEyeLeft = LeftEyeFrames.front();
			LeftEyeFrames.pop();
			LeftEyeFrames.push(eyePoses[ovrEye_Left]);
			eyePoses[ovrEye_Left] = currEyeLeft;

			currEyeRight = RightEyeFrames.front();
			RightEyeFrames.pop();
			RightEyeFrames.push(eyePoses[ovrEye_Right]);
			eyePoses[ovrEye_Right] = currEyeRight;
		}
	}
	

	// Rendering Delay
	if (Num_RenderingDelayFrames > 0) {
		if (NumRenderingFrames == Num_RenderingDelayFrames) {
			NumRenderingFrames = 0;
			currHandLeft = { handPoses[ovrHand_Left].Orientation,handPoses[ovrHand_Left].Position };
			currHandRight = { handPoses[ovrHand_Right].Orientation,handPoses[ovrHand_Right].Position };
			oldHandLeft = { currHandLeft.Orientation, currHandLeft.Position };
			oldHandRight = { currHandRight.Orientation, currHandRight.Position };

			currEyeLeft = { eyePoses[ovrEye_Left].Orientation, eyePoses[ovrEye_Left].Position };
			currEyeRight = { eyePoses[ovrEye_Right].Orientation, eyePoses[ovrEye_Right].Position };
			oldEyeLeft = { currEyeLeft.Orientation, currEyeLeft.Position };
			oldEyeRight = { currEyeRight.Orientation, currEyeRight.Position };
		}
		else {
			NumRenderingFrames++;
			
			handPoses[ovrHand_Left] = { oldHandLeft.Orientation, oldHandLeft.Position };
			handPoses[ovrHand_Right] = { oldHandRight.Orientation, oldHandRight.Position };
			eyePoses[ovrEye_Left] = { oldEyeLeft.Orientation, oldEyeLeft.Position };
			eyePoses[ovrEye_Right] = { oldEyeRight.Orientation, oldEyeRight.Position };

		}
	}
	else {
		currHandLeft = { handPoses[ovrHand_Left].Orientation,handPoses[ovrHand_Left].Position };
		currHandRight = { handPoses[ovrHand_Right].Orientation,handPoses[ovrHand_Right].Position };
		oldHandLeft = { currHandLeft.Orientation, currHandLeft.Position };
		oldHandRight = { currHandRight.Orientation, currHandRight.Position };

		currEyeLeft = { eyePoses[ovrEye_Left].Orientation, eyePoses[ovrEye_Left].Position };
		currEyeRight = { eyePoses[ovrEye_Right].Orientation, eyePoses[ovrEye_Right].Position };
		oldEyeLeft = { currEyeLeft.Orientation, currEyeLeft.Position };
		oldEyeRight = { currEyeRight.Orientation, currEyeRight.Position };
	}


	

    int curIndex;
    ovr_GetTextureSwapChainCurrentIndex(_session, _eyeTexture, &curIndex);
    GLuint curTexId;
    ovr_GetTextureSwapChainBufferGL(_session, _eyeTexture, curIndex, &curTexId);
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, _fbo);
    glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, curTexId, 0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    ovr::for_each_eye([&](ovrEyeType eye)
    {
      const auto& vp = _sceneLayer.Viewport[eye];
      glViewport(vp.Pos.x, vp.Pos.y, vp.Size.w, vp.Size.h);
      _sceneLayer.RenderPose[eye] = originEyePoses[eye];

	 
	  // HEADING TRACKING MODE SWITCH
	  if (headMode == REGULAR) {
		  currEyeLeft = { eyePoses[ovrEye_Left].Orientation, eyePoses[ovrEye_Left].Position };
		  currEyeRight = { eyePoses[ovrEye_Right].Orientation, eyePoses[ovrEye_Right].Position };
		  oldEyeLeft = { currEyeLeft.Orientation, currEyeLeft.Position };
		  oldEyeRight = { currEyeRight.Orientation, currEyeRight.Position };
	  }
	  else if (headMode == ORIENTATION_ONLY) {
		  currEyeLeft = { eyePoses[ovrEye_Left].Orientation, oldEyeLeft.Position };
		  currEyeRight = { eyePoses[ovrEye_Right].Orientation, oldEyeRight.Position };
	  }
	  else if (headMode == POSITION_ONLY) {
		  currEyeLeft = { oldEyeLeft.Orientation, eyePoses[ovrEye_Left].Position };
		  currEyeRight = { oldEyeRight.Orientation, eyePoses[ovrEye_Right].Position };
	  }
	  else if (headMode == NO_TRACKING) {
		  currEyeLeft = { oldEyeLeft.Orientation, oldEyeLeft.Position };
		  currEyeRight = { oldEyeRight.Orientation, oldEyeRight.Position };
	  }

	  
	  

	  // VIEWING MODE SWITCH
	  mat4 centerCameraMatrix;
	  switch (viewMode)
	  {
	  case STEREO:
		  if (eye == ovrEye_Left) {
			  renderScene(_eyeProjections[ovrEye_Left], ovr::toGlm(currEyeLeft), true);
		  }
		  else {
			  renderScene(_eyeProjections[ovrEye_Right], ovr::toGlm(currEyeRight), false);
		  }
		  break;
	  case MONO:
		  centerCameraMatrix = glm::scale((ovr::toGlm(currEyeLeft) + ovr::toGlm(currEyeRight)),glm::vec3(0.5f));
		  renderScene(_eyeProjections[eye], centerCameraMatrix, true);
		  break;
	  case LEFT_EYE_ONLY:
		  if (eye == ovrEye_Left) {
			  renderScene(_eyeProjections[ovrEye_Left], ovr::toGlm(currEyeLeft), true);
		  }
		  break;
	  case RIGHT_EYE_ONLY:
		  if (eye == ovrEye_Right) {
			  renderScene(_eyeProjections[ovrEye_Right], ovr::toGlm(currEyeRight), false);
		  }
		  break;
	  case INVERTED_EYES:
		  if (eye == ovrEye_Left) {
			  renderScene(_eyeProjections[ovrEye_Right], ovr::toGlm(currEyeRight), false);
		  }
		  else {
			  renderScene(_eyeProjections[ovrEye_Left], ovr::toGlm(currEyeLeft), true);
		  }
		  break;
	  }

	  // IOD
	  float horizontal_center = (_viewScaleDesc.HmdToEyePose[ovrEye_Left].Position.x + _viewScaleDesc.HmdToEyePose[ovrEye_Right].Position.x) / 2;
	
	  _viewScaleDesc.HmdToEyePose[ovrEye_Left].Position.x = horizontal_center - IOD / 2;
	  _viewScaleDesc.HmdToEyePose[ovrEye_Right].Position.x = horizontal_center + IOD / 2;

    });
    glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, 0, 0);
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
    ovr_CommitTextureSwapChain(_session, _eyeTexture);
    ovrLayerHeader* headerList = &_sceneLayer.Header;
    ovr_SubmitFrame(_session, frame, &_viewScaleDesc, &headerList, 1);

    GLuint mirrorTextureId;
    ovr_GetMirrorTextureBufferGL(_session, _mirrorTexture, &mirrorTextureId);
    glBindFramebuffer(GL_READ_FRAMEBUFFER, _mirrorFbo);
    glFramebufferTexture2D(GL_READ_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, mirrorTextureId, 0);
    glBlitFramebuffer(0, 0, _mirrorSize.x, _mirrorSize.y, 0, _mirrorSize.y, _mirrorSize.x, 0, GL_COLOR_BUFFER_BIT,
                      GL_NEAREST);
    glBindFramebuffer(GL_READ_FRAMEBUFFER, 0);
  }

  virtual void renderScene(const glm::mat4& projection, const glm::mat4& headPose, bool isLeft) = 0;
};

//////////////////////////////////////////////////////////////////////
//
// The remainder of this code is specific to the scene we want to 
// render.  I use glfw to render an array of cubes, but your 
// application would perform whatever rendering you want
//

#include <vector>
#include <algorithm>
#include "shader.h"
#include "Cube.h"
#include<stdio.h> 
#include "Model.h"
#include "Mesh.h"
#include "Line.h"


// a class for building and rendering cubes
class Scene
{

  // Grid Size
  const unsigned int GRID_SIZE{ 10 };

  // Program
  std::vector<glm::mat4> instance_positions;
  GLuint instanceCount;
  GLuint shaderID;
  GLuint lineShaderID;
  GLuint skyboxShaderID;
  GLuint highlightShaderID;

  //std::unique_ptr<TexturedCube> cube;
  std::unique_ptr<Skybox> skybox_righteye;
  std::unique_ptr<Skybox> skybox_lefteye;

  // Final Projects

  // Cubes
  // My Board
  std::unique_ptr<TexturedCube> my_board_cube_normal;
  std::unique_ptr<TexturedCube> my_board_cube_A;
  std::unique_ptr<TexturedCube> my_board_cube_B;
  std::unique_ptr<TexturedCube> my_board_cube_C;
  std::unique_ptr<TexturedCube> my_board_cube_S;
  std::unique_ptr<TexturedCube> my_board_cube_P;
  std::unique_ptr<TexturedCube> my_board_cube_occupied;
  std::unique_ptr<TexturedCube> my_board_cube_missed;
  std::unique_ptr<TexturedCube> my_board_cube_shooted;

  std::vector<glm::vec3> my_board_positions;
  std::vector<glm::vec3> rival_board_positions;
  std::vector<glm::mat4> my_board_matrices;
  std::vector<glm::mat4> rival_board_matrices;

  // Lines
  std::unique_ptr<Line> line;

  // Game Data
  int myBoard[10][10];
  int rivalBoard[10][10];
  

  //
  const unsigned int EMPTY = 0;
  const unsigned int OCCUPIED = 1;
  //
  const unsigned int TAKEN = 3;
  const unsigned int MISSED = 4;
  const unsigned int SHOOTED = 5;

 

  //

  // Extra Credit 1
  std::unique_ptr<Skybox> skybox_customized_1;
  std::unique_ptr<Skybox> skybox_customized_2;

  

public:
	//
	int x_cord, y_cord;

	// Warships
	std::vector<std::pair<int, int>> A; // Aircraft_carrier * 1 (5*1 grid)
	std::vector<std::pair<int, int>> B; // Battleship * 1 (4*1 grid)
	std::vector<std::pair<int, int>> C; // Cruiser * 1 (3*1 grid)
	std::vector<std::pair<int, int>> S; // Submarine * 1 (2*1 grid)
	std::vector<std::pair<int, int>> P; // Patrol * 1 (2*1 grid)
	std::vector<std::pair<int, int>> All; // Patrol * 1 (2*1 grid)

	// Cube Size: Default Size = 30cm. Support a range from 0.01m to 0.5m.
	float scale_ratio;
	const float MAX_CUBE_SIZE = 0.25f; // 50cm
	const float MIN_CUBE_SIZE = 0.005f; // 1cm
	const float DEFAULT_CUBE_SIZE = 0.15f; // 30cm

  Scene()
  {
    // Create two cube
    instance_positions.push_back(glm::translate(glm::mat4(1.0f), glm::vec3(0, 0, -0.3)));
    instance_positions.push_back(glm::translate(glm::mat4(1.0f), glm::vec3(0, 0, -0.9)));

    instanceCount = instance_positions.size();

    // Shader Program 
    skyboxShaderID = LoadShaders("skybox.vert", "skybox.frag");
	lineShaderID = LoadShaders("line.vert", "line.frag");
	highlightShaderID = LoadShaders("highlight.vert", "highlight.frag");

	shaderID = skyboxShaderID;

	// 10m wide sky box: size doesn't matter though
    skybox_righteye = std::make_unique<Skybox>("skybox_righteye");
	skybox_lefteye = std::make_unique<Skybox>("skybox_lefteye");

	skybox_righteye->toWorld = glm::scale(glm::mat4(1.0f), glm::vec3(5.0f));
	skybox_lefteye->toWorld = glm::scale(glm::mat4(1.0f), glm::vec3(5.0f));

	// Customized Skybox
	skybox_customized_1 = std::make_unique<Skybox>("skybox_customized_1");
	skybox_customized_1->toWorld = glm::scale(glm::mat4(1.0f), glm::vec3(5.0f));

	skybox_customized_2 = std::make_unique<Skybox>("skybox_customized_2");
	skybox_customized_2->toWorld = glm::scale(glm::mat4(1.0f), glm::vec3(5.0f));

	// Default Cube Size: Scale to 30cm: 0.15
	scale_ratio = DEFAULT_CUBE_SIZE;


	// Final Projects

	// Cube
	// My Board
	my_board_cube_normal = std::make_unique<TexturedCube>("My_Board_Normal_Square");
	my_board_cube_A = std::make_unique<TexturedCube>("My_Board_A_Square");
	my_board_cube_B = std::make_unique<TexturedCube>("My_Board_B_Square");
	my_board_cube_C = std::make_unique<TexturedCube>("My_Board_C_Square");
	my_board_cube_S = std::make_unique<TexturedCube>("My_Board_S_Square");
	my_board_cube_P = std::make_unique<TexturedCube>("My_Board_P_Square");
	my_board_cube_occupied = std::make_unique<TexturedCube>("My_Board_Occupied_Square");

	//10*10 Board Cube Positions 
	for (unsigned int z = 0; z < GRID_SIZE; ++z) {
		for (unsigned int x = 0; x < GRID_SIZE; ++x) {
			// My Board Position
			float xpos = -1.0f + 0.2f * x;
			float ypos = -1.0f;
			float zpos = -1.0f + 0.2f * z;
			glm::vec3 relativePosition = glm::vec3(xpos, ypos, zpos);
			my_board_positions.push_back(relativePosition);
			my_board_matrices.push_back(glm::translate(glm::mat4(1.0f), relativePosition));

			// Rival Board Position
			xpos = 0.04f * x;
			ypos = 0.0f;
			zpos = 0.04f * z;
			relativePosition = glm::vec3(xpos, ypos, zpos);
			rival_board_positions.push_back(relativePosition);
			rival_board_matrices.push_back(glm::translate(glm::mat4(1.0f), relativePosition));
		}
	}

	// Line
	line = std::make_unique<Line>();

	


  }

  void render(const glm::mat4& projection, const glm::mat4& view, bool isLeft)
  {
   

    // Render Skybox : remove view translation
	if (skyboxMode == SKYBOX_ENTIRE) {
		if (isLeft) {
			skybox_lefteye->draw(shaderID, projection, view);
		}
		else {
			skybox_righteye->draw(shaderID, projection, view);
		}

		for (int i = 0; i < instanceCount; i++)
		{
			//cube->toWorld = instance_positions[i] * glm::scale(glm::mat4(1.0f), glm::vec3(scale_ratio));
			//cube->draw(shaderID, projection, view);
		}
	}
	else if (skyboxMode == SKYBOX_STEREO) {
		if (isLeft) {
			skybox_lefteye->draw(shaderID, projection, view);
		}
		else {
			skybox_righteye->draw(shaderID, projection, view);
		}
	}
	else if (skyboxMode == SKYBOX_MONO) {
		skybox_lefteye->draw(shaderID, projection, view);
	}
	else if(skyboxMode == SKYBOX_MONO_CUSTOMIZED_1) {
		skybox_customized_1->draw(shaderID, projection, view);
	}
	else if (skyboxMode == SKYBOX_MONO_CUSTOMIZED_2) {
		skybox_customized_2->draw(shaderID, projection, view);
	}

	// Final Project
	int row;
	int column;
	std::pair<int, int> curr;
	// My Board
	for (int i = 0; i < GRID_SIZE*GRID_SIZE; i++)
	{
		row = i % GRID_SIZE;
		column = (i - row) / GRID_SIZE;
		curr = std::make_pair(row, column);

	   if (column == y_cord && row == x_cord) {
		   if (warshipMode == A_Ship) {
			   my_board_cube_A->toWorld = my_board_matrices[i] * glm::scale(glm::mat4(1.0f), glm::vec3(0.1f, 0.005f, 0.1f));
			   my_board_cube_A->draw(skyboxShaderID, projection, view);
		   }
		   else if (warshipMode == B_Ship) {
			   my_board_cube_B->toWorld = my_board_matrices[i] * glm::scale(glm::mat4(1.0f), glm::vec3(0.1f, 0.005f, 0.1f));
			   my_board_cube_B->draw(skyboxShaderID, projection, view);
		   }
		   else if (warshipMode == C_Ship) {
			   my_board_cube_C->toWorld = my_board_matrices[i] * glm::scale(glm::mat4(1.0f), glm::vec3(0.1f, 0.005f, 0.1f));
			   my_board_cube_C->draw(skyboxShaderID, projection, view);
		   }
		   else if (warshipMode == S_Ship) {
			   my_board_cube_S->toWorld = my_board_matrices[i] * glm::scale(glm::mat4(1.0f), glm::vec3(0.1f, 0.005f, 0.1f));
			   my_board_cube_S->draw(skyboxShaderID, projection, view);
		   }
		   else if (warshipMode == P_Ship) {
			   my_board_cube_P->toWorld = my_board_matrices[i] * glm::scale(glm::mat4(1.0f), glm::vec3(0.1f, 0.005f, 0.1f));
			   my_board_cube_P->draw(skyboxShaderID, projection, view);
		   }
		}
	   if (std::find(A.begin(), A.end(), curr) != A.end() || std::find(B.begin(), B.end(), curr) != B.end() ||
		   std::find(C.begin(), C.end(), curr) != C.end() || std::find(S.begin(), S.end(), curr) != S.end() ||
		   std::find(P.begin(), P.end(), curr) != P.end()) {
		   my_board_cube_occupied->toWorld = my_board_matrices[i] * glm::scale(glm::mat4(1.0f), glm::vec3(0.1f, 0.005f, 0.1f));
		   my_board_cube_occupied->draw(skyboxShaderID, projection, view);
	   }
		else {
			my_board_cube_normal->toWorld = my_board_matrices[i] * glm::scale(glm::mat4(1.0f), glm::vec3(0.1f, 0.005f, 0.1f));
			my_board_cube_normal->draw(skyboxShaderID, projection, view);
		}
	}

	// Rival Board
	glm::mat4 LHOrientationPosition;
	// Translation
	glm::mat4 translate = glm::translate(glm::mat4(1.0f), vec3(handPosition[ovrHand_Left].x, handPosition[ovrHand_Left].y, handPosition[ovrHand_Left].z));
	// Rotation
	float x = handRotation[ovrHand_Left].x;
	float y = handRotation[ovrHand_Left].y;
	float z = handRotation[ovrHand_Left].z;
	float w = handRotation[ovrHand_Left].w;
	// Quaternion to Rotation Matrix
	glm::mat4 rotate = glm::mat4(1 - 2 * y*y - 2 * z*z, 2 * x*y + 2 * w*z, 2 * x*z - 2 * w*y, 0.0f,
		2 * x*y - 2 * w*z, 1 - 2 * x*x - 2 * z*z, 2 * y*z + 2 * w*x, 0.0f,
		2 * x*z + 2 * w*y, 2 * y*z - 2 * w*x, 1 - 2 * x*x - 2 * y*y, 0.0f,
		0.0f, 0.0f, 0.0f, 1.0f);
	// Left-Hand orientation and position
	LHOrientationPosition = translate * rotate;


	//
	 
	for (int i = 0; i < GRID_SIZE*GRID_SIZE; i++)
	{
		my_board_cube_normal->toWorld = LHOrientationPosition * rival_board_matrices[i] * glm::scale(glm::mat4(1.0f), glm::vec3(0.02f, 0.005f, 0.02f));
		my_board_cube_normal->draw(shaderID, projection, view);
	}


    //
	/*float delta = (-1.0f - handPosition[ovrHand_Right].y) / handRotation[ovrHand_Right].y;
	glm::vec3 endPos = glm::vec3(handPosition[ovrHand_Right].x + delta * handRotation[ovrHand_Right].x, -1.0f, handPosition[ovrHand_Right].z + delta * handRotation[ovrHand_Right].z);
	float minDistSq = FLT_MAX;
	float currDistSq;
	int closetIdx;

	glUseProgram(lineShaderID);
	for (int i = 0; i < GRID_SIZE*GRID_SIZE; i++)
	{
		currDistSq = (my_board_positions[i].x - endPos.x)*(my_board_positions[i].x - endPos.x) + (my_board_positions[i].y - endPos.y)*(my_board_positions[i].y - endPos.y) +
			(my_board_positions[i].z - endPos.z)*(my_board_positions[i].z - endPos.z);
		if (currDistSq < minDistSq) {
			minDistSq = currDistSq;
			closetIdx = i;
		}
	}

	line->update(vec3(handPosition[ovrHand_Right].x, handPosition[ovrHand_Right].y, handPosition[ovrHand_Right].z), my_board_positions[closetIdx]);
	line->draw(lineShaderID, projection, view);*/
	
  }

  void placeWarship(std::vector<std::pair<int, int>> & currWarship, unsigned int size) {
	  switch (directionMode) {
	  case NONE:
		  return;
	  case UP:
		  if (y_cord >= (size-1)) {
			  currWarship.clear();
			  for (unsigned int i = 0; i < size; i++) {
				  if (checkAvailability(std::make_pair(x_cord, y_cord - i))) {
					  currWarship.push_back(std::make_pair(x_cord, y_cord - i));
				  }
				  else {
					  currWarship.clear();
					  break;
				  }
				  
			  } 
		  }
		  break;
	  case RIGHT:
		  if ((9 - x_cord) >= (size - 1)) {
			  currWarship.clear();
			  for (unsigned int i = 0; i < size; i++) {
				  if (checkAvailability(std::make_pair(x_cord + i, y_cord))) {
					  currWarship.push_back(std::make_pair(x_cord + i, y_cord));
				  }
				  else {
					  currWarship.clear();
					  break;
				  }
			  }
		  }
		  break;
	  case DOWN:
		  if ((9 - y_cord) >= (size - 1)) {
			  currWarship.clear();
			  for (unsigned int i = 0; i < size; i++) {
				  if (checkAvailability(std::make_pair(x_cord, y_cord + i))) {
					  currWarship.push_back(std::make_pair(x_cord, y_cord + i));
				  }
				  else {
					  currWarship.clear();
					  break;
				  }
			  }
		  }
		  break;
	  case LEFT:
		  if (x_cord >= (size - 1)) {
			  currWarship.clear();
			  for (unsigned int i = 0; i < size; i++) {
				  if (checkAvailability(std::make_pair(x_cord - i, y_cord))) {
					  currWarship.push_back(std::make_pair(x_cord - i, y_cord));
				  }
				  else {
					  currWarship.clear();
					  break;
				  }		  
			  }
		  }
		  break;
	  }
  }

  bool checkAvailability(std::pair<int,int> curr) {
	  if (std::find(A.begin(), A.end(), curr) != A.end() || std::find(B.begin(), B.end(), curr) != B.end() ||
		  std::find(C.begin(), C.end(), curr) != C.end() || std::find(S.begin(), S.end(), curr) != S.end() ||
		  std::find(P.begin(), P.end(), curr) != P.end()) {
		  return false;
	  }
	  return true;
  }

};

class Cursor {

	// Shader ID
	GLuint shaderID;

	// Cursor
	std::unique_ptr<Model> cursor;

	// User's Dominant Hand's Controller Position 
	vec3 position;


public:
	Cursor() {
		shaderID = LoadShaders("shader.vert", "shader.frag");
		cursor = std::make_unique<Model>("webtrcc.obj");
	}

	/* Render sphere at User's Dominant Hand's Controller Position */
	void render(const glm::mat4& projection, const glm::mat4& view, vec3 pos) {
		position = pos;
		glm::mat4 toWorld = glm::translate(glm::mat4(1.0f), position) * glm::scale(glm::mat4(1.0f), glm::vec3(0.02f));
		cursor->Draw(shaderID, projection, view, toWorld);
	}

};


// An example application that renders a simple cube
class ExampleApp : public RiftApp
{
  std::shared_ptr<Scene> scene;

  // Thumb Sticks
  bool thumbStickRightPushed_X;
  bool thumbStickRightPushed_Y;
  

  // Buttons Tracking
  bool buttonXPressed;
  bool buttonAPressed;
  bool buttonBPressed;
  bool buttonYPressed;

  // Cursor
  std::unique_ptr<Cursor> cursor;

  // Triggers Tracking
  bool triggerLeftIndexClicked;
  bool triggerRightIndexClicked;
  bool triggerLeftMiddleClicked;
  bool triggerRightMiddleClicked;

public:

	
  ExampleApp()
  {
  }

protected:
  void initGl() override
  {
    RiftApp::initGl();
    glClearColor(0.2f, 0.2f, 0.2f, 0.0f);
    glEnable(GL_DEPTH_TEST);
    ovr_RecenterTrackingOrigin(_session);
    scene = std::shared_ptr<Scene>(new Scene());

	// Modes
	gameMode = PREPARE;
	// Prepare Mode
	directionMode = NONE;
	warshipMode = A_Ship;

	// Thumb Stick
	thumbStickRightPushed_X = false;
	thumbStickRightPushed_Y = false;
	

	// Buttons
	buttonXPressed = false;
	buttonAPressed = false;
	buttonBPressed = false;
	buttonYPressed = false;

	// Modes
	headMode = REGULAR;
	viewMode = STEREO;
	skyboxMode = SKYBOX_ENTIRE;

	// Cursor
	cursor = std::unique_ptr<Cursor>(new Cursor());

	// Triggers
	triggerLeftIndexClicked = false;
	triggerRightIndexClicked = false;
	triggerLeftMiddleClicked = false;
	triggerRightMiddleClicked = false;

	scene->x_cord = 0;
	scene->y_cord = 0;

  }

  void shutdownGl() override
  {
    scene.reset();
  }

  void renderScene(const glm::mat4& projection, const glm::mat4& headPose, bool isLeft) override
  {
	  
	  //
	  if (gameMode == PREPARE) {
		  if (OVR_SUCCESS(ovr_GetInputState(_session, ovrControllerType_Touch, &inputState))) {


			  // The LEFT INDEX trigger reduces the X_Coordinate of Selected Grid
			  if (inputState.IndexTrigger[ovrHand_Left] > 0.5f)
			  {
				  if (!triggerLeftIndexClicked) {
					  triggerLeftIndexClicked = true;
					  if (scene->x_cord > 0) { scene->x_cord--; UpdateWarshipPosition();
					  }
				  }
			  }
			  else { triggerLeftIndexClicked = false; }

			  // The RIGHT INDEX trigger increases the X_Coordinate of Selected Grid
			  if (inputState.IndexTrigger[ovrHand_Right] > 0.5f)
			  {
				  if (!triggerRightIndexClicked) {
					  triggerRightIndexClicked = true;
					  if (scene->x_cord < 9) { scene->x_cord++; UpdateWarshipPosition();
					  }
				  }
			  }
			  else { triggerRightIndexClicked = false; }

			  // The LEFT MIDDLE finger trigger reduces the Y_Coordinate of Selected Grid
			  if (inputState.HandTrigger[ovrHand_Left] > 0.5f)
			  {
				  if (!triggerLeftMiddleClicked) {
					  triggerLeftMiddleClicked = true;
					  if (scene->y_cord > 0) { scene->y_cord--; UpdateWarshipPosition();
					  }
				  }
			  }
			  else { triggerLeftMiddleClicked = false; }

			  // The RIGHT MIDDLE finger trigger increases the Y_Coordinate of Selected Grid
			  if (inputState.HandTrigger[ovrHand_Right] > 0.5f)
			  {
				  if (!triggerRightMiddleClicked) {
					  triggerRightMiddleClicked = true;
					  if (scene->y_cord < 9) { scene->y_cord++; UpdateWarshipPosition();
					  }
				  }
			  }
			  else { triggerRightMiddleClicked = false; }

			  // Button B
			  if (inputState.Buttons & ovrButton_B)
			  {
				  if (!buttonBPressed) {
					  buttonBPressed = true;
					  UpdateButtonB();
				  }
			  }
			  else { buttonBPressed = false; }

			  // Button A
			  if (inputState.Buttons & ovrButton_A)
			  {
				  if (!buttonAPressed) {
					  buttonAPressed = true;
					  UpdateButtonA();
				  }
			  }
			  else { buttonAPressed = false; }

			  // Button X
			  if (inputState.Buttons & ovrButton_X)
			  {
				  if (!buttonXPressed) {
					  buttonXPressed = true;
					  UpdateButtonX();
				  }
			  }
			  else { buttonXPressed = false; }
	      } 
		 
	  }

	  
	  
	  // Vary the physical size of BOTH CUBES with the left thumb stick left/right, changing the size of the cubes from 30cm to smaller or larger values. 
	  if (OVR_SUCCESS(ovr_GetInputState(_session, ovrControllerType_Touch, &inputState)))
	  {
		  if (inputState.Thumbstick[ovrHand_Left].x > 0.0f) // MOVE RIGHT
		  {
			  if (scene->scale_ratio > scene->MIN_CUBE_SIZE) {
				  scene->scale_ratio -= 0.0001f;
			  }
		  }
		  else if (inputState.Thumbstick[ovrHand_Left].x < 0.0f) { // MOVE LEFT
			  if (scene->scale_ratio < scene->MAX_CUBE_SIZE) {
				  scene->scale_ratio += 0.0001f;
			  }
		  }
	  }
	  // PUSH down on the thumb stick to reset the cubes to their initial sizes (30cm)
	  if (OVR_SUCCESS(ovr_GetInputState(_session, ovrControllerType_Touch, &inputState))) {
		  if (inputState.Buttons & ovrButton_LThumb)
		  {
			  scene->scale_ratio = scene->DEFAULT_CUBE_SIZE;
		  }
	  }

	

	 



	  
	  // Extra Credit 2: Smooth Controller (Button Y)
	  if (OVR_SUCCESS(ovr_GetInputState(_session, ovrControllerType_Touch, &inputState))) {
		  if (inputState.Buttons & ovrButton_Y)
		  {
			  if (!buttonYPressed) {
				  buttonYPressed = true;
				  UpdateButtonY();
			  }
		  }
		  else {
			  buttonYPressed = false;
		  }
	  }

	  // Render a sphere at your dominant hand's controller position
	  cursor->render(projection, glm::inverse(headPose), vec3(currHandRight.Position.x, currHandRight.Position.y, currHandRight.Position.z));

	  // Render Scene
      scene->render(projection, glm::inverse(headPose), isLeft);

	
  }

  

  void UpdateButtonX() {
	  if (selectingMode == DONE) { selectingMode == SELECTING; std::cout << "Slecting" << std::endl; }
	  else if (selectingMode == SELECTING) { selectingMode == DONE; std::cout << "Done" << std::endl; }
  }

  void UpdateButtonA() {
	  if (warshipMode == A_Ship) { warshipMode = B_Ship; std::cout << "B" << std::endl; }
	  else if (warshipMode == B_Ship) { warshipMode = C_Ship; std::cout << "C" << std::endl; }
	  else if (warshipMode == C_Ship) { warshipMode = S_Ship; std::cout << "S" << std::endl; }
	  else if (warshipMode == S_Ship) { warshipMode = P_Ship; std::cout << "P" << std::endl; }
	  else if (warshipMode == P_Ship) { warshipMode = A_Ship; std::cout << "A" << std::endl; }
  }

  
  void UpdateButtonB() {
	  if (directionMode == NONE) { directionMode = UP; std::cout << "UP" << std::endl;}
	  else if (directionMode == UP) { directionMode = RIGHT; std::cout << "RIGHT" << std::endl;}
	  else if (directionMode == RIGHT) { directionMode = DOWN; std::cout << "DOWN" << std::endl;}
	  else if (directionMode == DOWN) { directionMode = LEFT; std::cout << "LEFT" << std::endl;}
	  else if (directionMode == LEFT) { directionMode = UP; std::cout << "UP" << std::endl;}

	  UpdateWarshipPosition();
	 // std::cout << scene->A.size() << std::endl;
	  //for (int i = 0; i < scene->A.size(); i++) {
		 // std::pair<int, int> ele = scene->A.at(i);
		  //std::cout << ele.first << " " << ele.second << std::endl;
	  //}
  }



  void UpdateButtonY() {
	  if (NumAverageMovingFrames < 45) {
		  NumAverageMovingFrames++;
		  std::cout << "Smoothing Mode : " << NumAverageMovingFrames << " frames" << std::endl;
	  }
  }

  void UpdateWarshipPosition() {
	  switch (warshipMode) {
	  case A_Ship:
		  scene->placeWarship(scene->A, 5);
		  break;
	  case B_Ship:
		  scene->placeWarship(scene->B, 4);
		  break;
	  case C_Ship:
		  scene->placeWarship(scene->C, 3);
		  break;
	  case S_Ship:
		  scene->placeWarship(scene->S, 2);
		  break;
	  case P_Ship:
		  scene->placeWarship(scene->P, 2);
		  break;
	  }

  

  }
};

// Execute our example class
int main(int argc, char** argv)
{
  int result = -1;

  if (!OVR_SUCCESS(ovr_Initialize(nullptr)))
  {
    FAIL("Failed to initialize the Oculus SDK");
  }
  result = ExampleApp().run();

  ovr_Shutdown();
  return result;
}
