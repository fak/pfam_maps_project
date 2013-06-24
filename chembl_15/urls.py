from django.conf.urls.defaults import patterns, include, url

# Uncomment the next two lines to enable the admin:
from django.contrib import admin
admin.autodiscover()

urlpatterns = patterns('chembl_15.views',
    # Examples:
    url(r'^$', 'index', name='index'),
    url(r'^evidence/$', 'evidence_portal', name='evidence_portal'),
    url(r'^evidence/(?P<pfam_name>[\w-]+)/$', 'evidence', name='evidence'),
    url(r'^conflicts/$', 'conflicts_portal', name='conflicts_portal'),
    url(r'^conflicts/(?P<conflict_id>.+)/$', 'conflicts', name = 'conflicts'),
    url(r'^vote/conflicts/(?P<act>\d+)/(?P<conflict_id>.+)/$', 'vote' , name = 'vote'),
    url(r'^details/conflicts/(?P<act>\d+)/$', 'details', name = 'details'),
    )

urlpatterns += patterns('',
            url(r'^admin/', include(admin.site.urls)),
            )

