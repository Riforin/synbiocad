from django.conf.urls import url

from . import views

urlpatterns = [
	url(r'^$', views.index, name = 'index'),
	url(r'^charlocus/$', views.charLocus, name = 'charloc'),
	url(r'^charlocus/prebuilt/$', views.charPreBuilt, name = 'charpre'),
	url(r'^charlocus/standard/$', views.charStdORF, name = 'charstd'),
	url(r'^existinglocus/$', views.existingLocus, name = 'charexist'),
	url(r'^existinglocus/delete/$', views.existingDelete, name = 'existdel'),
	url(r'^existinglocus/prebuilt/$', views.existingPremade, name = 'existpre'),
	url(r'^existinglocus/standard/$', views.existingStd, name = 'existstd'),
	url(r'^existinglocus/custom/$', views.existingCustom, name = 'existcustom'),
	url(r'^existinglocus/nearby/$', views.existingNearby, name = 'existnear'),
	url(r'^customlocus/$', views.customLocus, name = 'customloc'),
	url(r'^customlocus/nonvary/$', views.customNonVary, name = 'customnonvary'),
	url(r'^customlocus/vary/$', views.customVary, name = 'customvary'),
	url(r'^results/$', views.results, name = 'results')
]