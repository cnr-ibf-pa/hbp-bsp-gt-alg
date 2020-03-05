from django.urls import path
from django.conf.urls import url, include
from django.conf import settings
from django.conf.urls.static import static
from django.contrib.auth import views as auth_views
from django.contrib import admin

from . import views

urlpatterns = [
	url(r'^$', views.index),    
    url(r'^actionpage/$', views.actionpage),
    url(r'^log_user/(?P<user_id>[0-9]+)/(?P<pagename>[0-9a-zA-Z_\-]+)$', views.log_user),
]
