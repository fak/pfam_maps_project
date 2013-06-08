from django.conf.urls.defaults import patterns, include, url

# Uncomment the next two lines to enable the admin:
from django.contrib import admin
admin.autodiscover()

urlpatterns = patterns('chembl_15.views',
    # Examples:
    url(r'^$', 'index', name='index'),
    url(r'^evidence/(?P<pfam_name>[\w-]+)/$', 'evidence', name='evidence')
)
urlpatterns += patterns('',
            url(r'^admin/', include(admin.site.urls)),
            )

